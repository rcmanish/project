# BP35 PPG pulse-rate variability results

## Raw metrics

| phase | SDNN (ms) | RMSSD (ms) | LF (ms²) | HF (ms²) |
| --- | --- | --- | --- | --- |
| Pre | 69.32 | 64.88 | 1931.29 | 1425.40 |
| During | 76.71 | 59.34 | 1571.60 | 1330.85 |
| Post | 72.79 | 61.53 | 2640.58 | 1479.36 |

## Percentage changes

| index | Pre → During (%) | During → Post (%) |
| --- | --- | --- |
| SDNN (ms) | +10.65 | -5.10 |
| RMSSD (ms) | -8.53 | +3.69 |
| LF (ms²) | -18.62 | +68.02 |
| HF (ms²) | -6.63 | +11.16 |

## During-intervention segmentation

| Nominal activity | Start (s) | End (s) | Intervals |
| --- | --- | --- | --- |
| Listening 1 | 0 | 180 | 248 |
| Task while listening | 180 | 420 | 335 |
| Listening 2 | 420 | 600 | 254 |

## Data quality and interpretation

| phase | recorded_span_s | analyzed_span_s | detected_peaks | corrected_intervals | artifact_percent |
| --- | --- | --- | --- | --- | --- |
| Pre | 298.88 | 298.88 | 419 | 8 | 1.91 |
| During | 706.69 | 600.00 | 838 | 18 | 2.15 |
| Post | 298.40 | 298.40 | 412 | 7 | 1.70 |

* All three archives report device MD-V5-0000471, firmware 1.14.0, and a nominal PPG sampling rate of 25 Hz.

* Pre and Post recordings are shorter than 300 s; no samples were fabricated.
* During uses the first 600 s of a 706.7 s recording. The nominal task boundaries at 180 and 420 s are not verified by event markers.
* Metrics are PPG pulse-rate variability estimates, not ECG-derived HRV. Interpret comparisons cautiously, especially because During is twice as long as Pre and Post.
* The 25 Hz PPG sampling grid quantizes peak times in approximately 40 ms steps, limiting interval precision, especially for RMSSD.
* Artifact correction uses a 300–1500 ms physiologic range and a 20% centered 11-interval local-median threshold, followed by cubic interpolation.
* PSD uses 4 Hz cubic-spline interpolation and approximately 300 s Hann-window Welch segments with 50% overlap; LF and HF use exact-edge trapezoidal integration.

## Method references

* [Kubios HRV Scientific User's Guide](https://www.kubios.com/downloads/HRV-Scientific-Users-Guide.pdf): standard 4 Hz interpolation, LF/HF bands, and Welch spectral settings.
* [Task Force HRV standards (1996)](https://www.escardio.org/static-file/Escardio/Guidelines/Scientific-Statements/guidelines-Heart-Rate-Variability-FT-1996.pdf): conventional definitions of short-term HRV bands and metrics.
