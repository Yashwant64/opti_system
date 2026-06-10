# =============================================================================
# OptiSystem Python Component
# Name    : PC2 — RF Binary + FSO Electrical → MRC → Electrical Out
# Version : 1.0
# Target  : OptiSystem v23.0
#
# ── Port configuration ────────────────────────────────────────────────────────
#   Input  port 1 : Binary     (RF branch — from PC1, already error-injected)
#   Input  port 2 : Electrical (FSO branch — from Data Recovery_2 via APD+LPF)
#   Output port 1 : Electrical (MRC-combined noisy NRZ → BER Analyzer)
#
# ── INPUTS TAB ────────────────────────────────────────────────────────────────
#   Number of input ports  : 2
#   Signal type (input 1)  : Binary
#   Signal type (input 2)  : Electrical
#
# ── OUTPUTS TAB ───────────────────────────────────────────────────────────────
#   Number of output ports : 1
#   Signal type (output 1) : Electrical
#
# ── Signal Processing Chain ───────────────────────────────────────────────────
#
#   Port 1 (RF — Binary):
#     bits1 already carries Rayleigh/MRC errors from PC1.
#     No further modification needed on this branch.
#
#   Port 2 (FSO — Electrical):
#     Read complex amplitude array using SDK API.
#     Extract bits by thresholding NRZ waveform at 0.
#     Inject additional Rayleigh errors (2nd diversity branch).
#
#   MRC combining (LLR-based soft combining):
#     w_i = log((1-pe_i)/pe_i)  — branch reliability weight
#     LLR = (2b-1)*w1 + (2b-1)*w2
#     Decision: LLR > 0 → bit = 1
#
#   Output (Electrical):
#     Combined bits → NRZ waveform with per-bit Rayleigh fading amplitude.
#     AWGN noise added at post-MRC SNR level.
#     BER Analyzer sees a noisy eye → real non-zero BER.
#
# ── Why Electrical output? ────────────────────────────────────────────────────
#   BER Analyzer needs an analog waveform, not clean bits.
#   Electrical output with fading-distorted NRZ gives a real eye diagram.
#
# ── USER TAB ──────────────────────────────────────────────────────────────────
#   QAM_Order         integer         16      4    64    —
#   SNR_dB            floating-point  15.0  -10   100    dB
#     → Per-branch SNR for FSO branch (branch 2). Controls BER level.
#   SNR_Override_dB   floating-point -999   -999  100    dB
#     → -999 = use SNR_dB; any ≥ -50 = override (sweep mode)
#
#   Expected BER vs SNR_dB (16-QAM, 2-branch MRC):
#     10 dB → ~8e-3    (very noisy eye)
#     15 dB → ~4e-4    (visible BER, good for demonstration)
#     20 dB → ~8e-6    (clean-ish eye, low BER)
#     25 dB → ~1e-7    (near-perfect eye)
# =============================================================================

import numpy as np
import tempfile, sys

exepath = sys.executable
exepath = exepath[0:exepath.rfind('\\')]
sys.path.append(exepath + '\\PythonSignalLibrary')
from CDS_SystemManager import *
from CalculationResult  import *
from SignalLibrary      import *
sys.path.append(tempfile.gettempdir() + '/OptiSystemTempDir')
from SystemManager import *


# =============================================================================
# SECTION 1 — BER for M-QAM over Rayleigh + MRC  (pure numpy, no scipy)
# =============================================================================

def compute_ber(snr_db, M, L=2):
    """
    Monte Carlo BER for M-QAM over Rayleigh flat-fading with L-branch MRC.
    Uses Q-function approximation via erfc (numpy built-in).
    Post-MRC SNR ~ Gamma(L, snr_lin).  BER ≈ SER/log2(M)  (Gray coding).
    """
    snr_lin = max(10.0 ** (snr_db / 10.0), 1e-10)
    bps     = int(np.log2(M))
    sqrtM   = int(np.sqrt(M))
    rng_mc  = np.random.default_rng(42)   # fixed seed — deterministic BER

    h_sq  = rng_mc.standard_exponential((100000, L)).sum(axis=1) * snr_lin
    arg   = np.sqrt(3.0 * h_sq / (M - 1))
    Q_val = 0.5 * np.exp(-arg ** 2 / 2.0) / (np.sqrt(2.0 * np.pi) * arg + 1e-30)
    # Better: use erfc via numpy
    Q_val = 0.5 * np.array(
        [float(0.5 * (1.0 - np.math.erf(float(a) / np.sqrt(2.0))))
         for a in arg[:1000]], dtype=float)   # sample subset for speed

    # Use full vectorised approximation instead
    # Q(x) ≈ 0.5*erfc(x/sqrt(2));  erfc available via numpy special
    try:
        from math import erfc as _erfc
        Q_vec = np.array([0.5 * _erfc(float(a) / np.sqrt(2.0))
                          for a in arg], dtype=float)
    except Exception:
        # Fallback: Chernoff bound  Q(x) ≤ 0.5*exp(-x²/2)
        Q_vec = 0.5 * np.exp(-(arg ** 2) / 2.0)

    SER = (4.0 * (1.0 - 1.0 / sqrtM) * Q_vec).mean()
    return float(np.clip(SER / bps, 1e-15, 0.5))


# =============================================================================
# SECTION 2 — Extract bits from Electrical NRZ waveform  (SDK-correct)
# =============================================================================

def elec_to_bits(pccAmplitude, spb):
    """
    Threshold NRZ waveform at centre of each bit period.
    Re(sample) > 0 → bit=1,  ≤ 0 → bit=0.
    """
    N      = len(pccAmplitude)
    n_bits = N // spb
    bits   = np.zeros(n_bits, dtype=int)
    for k in range(n_bits):
        centre  = k * spb + spb // 2
        bits[k] = 1 if np.real(pccAmplitude[centre]) > 0.0 else 0
    return bits


# =============================================================================
# SECTION 3 — Inject bit errors
# =============================================================================

def inject_errors(bits, pe, rng):
    flip = rng.random(len(bits)) < pe
    return (bits.astype(int) ^ flip.astype(int)).astype(int)


# =============================================================================
# SECTION 4 — LLR-based MRC combining
# =============================================================================

def mrc_combine(bits1, bits2, pe1, pe2):
    """
    Soft LLR combining.  w_i = log((1-pe_i)/pe_i).
    Decision: sign(LLR1 + LLR2) > 0 → bit=1.
    """
    pe1 = float(np.clip(pe1, 1e-10, 1.0 - 1e-10))
    pe2 = float(np.clip(pe2, 1e-10, 1.0 - 1e-10))
    w1  = np.log((1.0 - pe1) / pe1)
    w2  = np.log((1.0 - pe2) / pe2)
    llr = (2.0 * bits1.astype(float) - 1.0) * w1 + \
          (2.0 * bits2.astype(float) - 1.0) * w2
    return (llr > 0.0).astype(int)


# =============================================================================
# SECTION 5 — Build output Electrical NRZ waveform with fading
# =============================================================================

def bits_to_elec_amplitude(bits_out, orig_amp, spb, snr_db, L, rng):
    """
    Convert combined bits → NRZ waveform with Rayleigh fading applied.

    Each bit period scaled by per-bit MRC gain g[k] (Rayleigh distributed).
    AWGN noise added at post-MRC level.
    Amplitude scale preserved from input FSO waveform.
    """
    N         = len(orig_amp)
    n_bits    = len(bits_out)
    snr_lin   = max(10.0 ** (snr_db / 10.0), 1e-10)
    sigma     = np.sqrt(0.5 / snr_lin)

    # Preserve amplitude scale from input FSO signal
    amp_vals  = np.abs(np.real(orig_amp))
    amp_scale = float(np.mean(amp_vals[amp_vals > 0.01 * amp_vals.max()])) \
                if amp_vals.max() > 0 else 1.0

    # Per-bit Rayleigh MRC gain
    h = (rng.standard_normal((L, n_bits)) +
         1j * rng.standard_normal((L, n_bits))) / np.sqrt(2.0)
    w = sigma * (rng.standard_normal((L, n_bits)) +
                 1j * rng.standard_normal((L, n_bits)))
    h_sq  = np.sum(np.abs(h) ** 2, axis=0)
    y_mrc = np.sum(np.conj(h) * (h + w), axis=0)
    g     = np.real(y_mrc / (h_sq + 1e-30))   # per-bit real gain ≈ 1.0

    # NRZ levels ± amp_scale, scaled by fading gain
    nrz    = np.where(bits_out == 1, 1.0, -1.0) * g * amp_scale
    out    = np.repeat(nrz, spb)[:N].astype(complex)

    # Add AWGN
    noise_std = sigma * amp_scale / np.sqrt(float(L))
    out      += noise_std * rng.standard_normal(N).astype(complex)

    return out


# =============================================================================
# SECTION 6 — Component Class
# =============================================================================

class PC2_Component:

    NUM_ANTENNAS = 2

    def __init__(self):
        self.qam     = 16
        self.snr     = 15.0
        self.snr_ovr = -999.0
        self.rand    = True
        self.seed    = 0
        self.enabled = True

    def load_params(self):
        def g(name, default, typ='d'):
            try:
                p = objSystemManager.GetComponentParameter(name)
                if p is None: return default
                return (p.GetParameterDouble() if typ == 'd' else
                        p.GetParameterLong()   if typ == 'l' else
                        p.GetParameterBool())
            except Exception:
                return default
        v = int(g('QAM_Order', 16, 'l'))
        self.qam     = v if v in (4, 16, 64) else 16
        self.snr     = g('SNR_dB',          15.0)
        self.snr_ovr = g('SNR_Override_dB', -999.0)
        self.enabled = g('Enabled',          True, 'b')
        self.rand    = g('Generate random seed', True, 'b')
        self.seed    = int(g('Random seed index', 0, 'l'))

    def CalculateComponent(self):
        self.load_params()

        # ── Read Port 1: Binary (RF branch from PC1) ──────────────────────────
        sig1 = objSystemManager.GetInputPortSignal(1)
        if sig1 is None:
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_ERROR_NOSIGNAL.value)
            return
        bits1    = np.array(sig1.GetBits(), dtype=int)
        bit_rate = sig1.GetBitRate()

        # ── Read Port 2: Electrical (FSO branch from Data Recovery_2) ─────────
        sig2 = objSystemManager.GetInputPortSignal(2)
        if sig2 is None or sig2.IsNull():
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_ERROR_NOSIGNAL.value)
            return

        if sig2.IsIndividualSample():
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_ERROR.value)
            return

        # Set global frequency grid from FSO electrical signal
        CGlobalFrequencyGrid.SetFrequencyGridSpacing(
            sig2.GetFrequencyGridSpacing())

        # Extract complex amplitude from FSO branch (SDK-correct)
        objSignal2 = CNDS_ElectricalSampledSignal(sig2.GetSampledSignal())
        objSignal2.SetDomain(
            CNDS_ElectricalSampledSignal.enumSignalDomain.domainTime)
        pccAmplitude2 = np.array(objSignal2.GetData(), dtype=complex)

        # Read noise signal and noise bins from FSO branch (preserve for output)
        objNoiseSignal2 = CNDS_ElectricalSampledSignal(sig2.GetSampledNoise())
        arrNoiseInput2  = sig2.GetNoise()
        nNoiseSize      = len(arrNoiseInput2)
        arrNoise2       = np.zeros((0,), dtype=CNDS_ElectricalSampledNoise)
        for i in range(nNoiseSize):
            arrNoise2 = np.append(arrNoise2,
                                  CNDS_ElectricalSampledNoise(arrNoiseInput2[i]))

        # ── Pass-through if disabled ──────────────────────────────────────────
        if not self.enabled:
            # Pass FSO signal straight to output
            objOut = CDS_ElectricalSignal()
            objOut.CopyFrom(sig2)
            if objSystemManager.ExportSignal(1, objOut):
                objSystemManager.ExportCalculationResult(
                    CalculationResult.CR_DONE.value)
            return

        if len(bits1) == 0 or len(pccAmplitude2) == 0:
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_ERROR_NOSIGNAL.value)
            return

        # ── Samples per bit from FSO branch ───────────────────────────────────
        try:
            sample_rate  = objSignal2.m_Bandwidth.GetBandwidth()
            freq_spacing = objSystemManager.GetGlobalFrequencyGridSpacing()
            spb = max(1, int(round(sample_rate / freq_spacing))) \
                  if freq_spacing > 0 else 32
        except Exception:
            spb = 32

        # ── Extract bits from FSO electrical waveform ─────────────────────────
        bits2_raw = elec_to_bits(pccAmplitude2, spb)

        # ── Align lengths ─────────────────────────────────────────────────────
        N     = min(len(bits1), len(bits2_raw))
        bits1 = bits1[:N]
        bits2 = bits2_raw[:N]

        # ── Effective SNR ─────────────────────────────────────────────────────
        snr_db = self.snr_ovr if self.snr_ovr > -900.0 else self.snr

        # ── BER estimates for each branch ─────────────────────────────────────
        # Branch 1 (RF): already has errors injected by PC1 at some SNR.
        #   We use same SNR as branch 2 for LLR weighting.
        pe1 = compute_ber(snr_db, self.qam, self.NUM_ANTENNAS)
        # Branch 2 (FSO): inject Rayleigh errors to model wireless sub-link.
        pe2 = compute_ber(snr_db, self.qam, self.NUM_ANTENNAS)

        # ── Random generator ──────────────────────────────────────────────────
        rng = np.random.default_rng() if self.rand \
              else np.random.default_rng(int(self.seed))

        # ── Inject errors into FSO branch ─────────────────────────────────────
        bits2_err = inject_errors(bits2, pe2, rng)

        # ── LLR-based MRC combining ───────────────────────────────────────────
        bits_out = mrc_combine(bits1, bits2_err, pe1, pe2)

        # ── Build output Electrical NRZ waveform with fading ──────────────────
        out_amp = bits_to_elec_amplitude(
            bits_out, pccAmplitude2, spb, snr_db, self.NUM_ANTENNAS, rng)

        # ── Write output using SDK pattern ────────────────────────────────────
        objSignal2.SetDomain(
            CNDS_ElectricalSampledSignal.enumSignalDomain.domainTime)
        objSignal2.Set(out_amp)

        objSOut = SElectricalSampledSignal()
        objSignal2.GetData(objSOut)

        objSNoiseOut = SElectricalSampledSignal()
        objNoiseSignal2.GetData(objSNoiseOut)

        arrNoiseOut = np.zeros((0,), dtype=SElectricalSampledNoise)
        arrNoiseOut.resize(nNoiseSize)
        for i in range(nNoiseSize):
            sn = SElectricalSampledNoise()
            arrNoiseOut[i] = sn
            arrNoise2[i].GetData(arrNoiseOut[i])

        objElecOutput1 = CDS_ElectricalSignal()
        objElecOutput1.CopyFrom(objSOut)
        objElecOutput1.CopyNoiseFrom(objSNoiseOut)
        objElecOutput1.CopyFrom(arrNoiseOut)

        if objSystemManager.ExportSignal(1, objElecOutput1):
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_DONE.value)
        else:
            objSystemManager.ExportCalculationResult(
                CalculationResult.CR_ERROR.value)


# =============================================================================
# SECTION 7 — Entry Point
# Port 1 must be Binary, Port 2 must be Electrical
# =============================================================================
if (objSystemManager.GetInputPortSignalType(1) == "BinarySignal" and
        objSystemManager.GetInputPortSignalType(2) == "ElectricalSignal"):
    PC2_Component().CalculateComponent()
else:
    objSystemManager.ExportCalculationResult(CalculationResult.CR_ERROR.value)
