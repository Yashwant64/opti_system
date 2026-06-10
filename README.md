# Performance Evaluation of Hybrid (FSO-RF) Link with Maximum Ratio Combining (MRC)

## Overview

This project presents the design and performance evaluation of a Hybrid Free Space Optical (FSO) and Radio Frequency (RF) communication system using Maximum Ratio Combining (MRC). The objective is to improve communication reliability under varying atmospheric turbulence conditions by combining the advantages of both FSO and RF links.

The system is simulated using OptiSystem, and its performance is analyzed under weak, moderate, and strong atmospheric turbulence conditions.

## Objectives

* Design a Hybrid FSO-RF communication system.
* Implement Maximum Ratio Combining (MRC) at the receiver.
* Analyze system performance under different turbulence conditions.
* Evaluate Quality Factor (Q-Factor), Bit Error Rate (BER), and Eye Diagrams.
* Compare system performance for different transmission distances and attenuation levels.

## System Components

### FSO Link

* CW Laser Source
* NRZ Pulse Generator
* Optical Modulator
* FSO Channel
* APD Photodetector
* Low Pass Filter

### RF Link

* RF Transmitter
* RF Channel
* RF Receiver

### Diversity Combining

* Maximum Ratio Combining (MRC)

## Methodology

1. Generate random binary data.
2. Transmit data simultaneously through FSO and RF links.
3. Model atmospheric turbulence effects on the FSO channel.
4. Apply MRC at the receiver to combine received signals.
5. Measure BER and Q-Factor using BER Analyzer.
6. Compare performance under different turbulence scenarios.

## Performance Metrics

### Bit Error Rate (BER)

BER is used to evaluate the reliability of the communication system.

### Q-Factor

Q-Factor indicates signal quality and system performance.

### Eye Diagram

Eye diagrams are analyzed to observe signal distortion and noise effects.

## Simulation Scenarios

### Weak Turbulence

* Highest Q-Factor
* Lowest BER
* Wide eye opening

### Moderate Turbulence

* Reduced Q-Factor
* Increased BER
* Partial eye closure

### Strong Turbulence

* Lowest Q-Factor
* Highest BER
* Significant eye distortion

## Results

The simulation results demonstrate that:

* Weak turbulence provides the best communication performance.
* Strong turbulence significantly degrades signal quality.
* BER increases with increasing turbulence severity.
* Q-Factor decreases with increasing distance and attenuation.
* MRC improves overall system reliability by exploiting diversity gain from both FSO and RF links.

## Tools Used

* OptiSystem
* Python
* MATLAB (optional for post-processing)
* LaTeX (for thesis documentation)

## Repository Structure

```text
├── OptiSystem_Project/
├── Results/
│   ├── Weak_Turbulence/
│   ├── Moderate_Turbulence/
│   └── Strong_Turbulence/
├── Graphs/
├── Thesis/
├── Python_Scripts/
└── README.md
```

## Future Work

* Adaptive modulation techniques
* Machine learning-based channel estimation
* Advanced diversity combining schemes
* Performance analysis under different weather conditions (fog, rain, haze)

## Author

**Yashwant Sahu**
M.Tech (Communication Systems)
National Institute of Technology (NIT) Bhopal

## License

This project is intended for academic and research purposes.
