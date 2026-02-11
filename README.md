# Alamouti QPSK TX/RX (PN Preamble, CFO Compensation) — MATLAB

This repository contains a MATLAB-based transmitter and receiver chain for a QPSK link using an Alamouti (2×1) space–time block code (STBC). The workflow supports offline generation of two transmit IQ streams (Tx1/Tx2) exported as CSV files, and offline decoding from recorded IQ text/CSV files using PN-preamble synchronization and CFO estimation/compensation.

## Contents
- **Transmitter (TX):** PN preamble + payload (16-bit length header + ASCII message) → QPSK (π/4) → two-slot preamble (Tx1-only then Tx2-only) → Alamouti mapping → RRC pulse shaping → export `leano1.csv` and `leano2.csv`.
- **Receiver (RX):** robust IQ file parsing → DC removal + normalization → RRC matched filter → symbol-phase scan → PN correlation for frame start → CFO estimation from preamble → CFO compensation → slot-A/slot-B preamble inspection → channel estimation → Alamouti combining → QPSK demod (4 rotations) → decode LEN+MSG.
  
---

## Equipment
This project was primarily evaluated using **offline IQ recordings** (TXT/CSV).  
If RF playback/capture is used, typical setups include:
- AnaPico APVSG40-4 Vector Signal Generatoror IQ playback source (optional) for transmission
- Anritsu MS2720xA-0743  for capturing IQ samples (optional)
- PC/laptop running MATLAB 2025.b for processing



---

## Tools and Dependencies
- **MATLAB** (tested with Communication System Toolbox)
- **Communication System Toolbox** (for RRC filters and QPSK mod/demod objects)
- Operating system: Windows/Linux/macOS (any MATLAB-supported OS)

---

## Techniques Used (High-level)
- Root-Raised-Cosine (RRC) pulse shaping and matched filtering
- PN-sequence based preamble detection via **normalized correlation**
- **Carrier Frequency Offset (CFO)** estimation from preamble phase slope and symbol-domain compensation
- Alamouti (2×1) STBC mapping and linear combining
- QPSK demodulation with **π/2 rotation search** to resolve phase ambiguity
- Structured payload parsing: **16-bit length header + ASCII message**

##Key TX parameters (default in code):
- sampleRate = 25e6
- symbolRate = 2.5e6
- samplesPerSymbol = 10
- rolloff = 0.35
- filterSpan = 6
- QPSK phase offset: pi/4


