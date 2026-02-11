# CEST Phantom Analysis Pipeline - User Guide

## Overview

This repository now includes two comprehensive CEST analysis scripts with automatic tube detection:

1. **`phantom_CEST_analysis_B0corr.m`** - CEST analysis with B0 correction only
2. **`phantom_CEST_analysis_FULL.m`** - Complete analysis with T1/T2/B1/B0 corrections

Both scripts automatically detect 24 tubes in your phantom using the same algorithm from `phantom_t1.m` and `phantom_T2.m`.

## 📋 Your 24-Tube Phantom Configuration

Your phantom contains the following tubes (labels will be applied automatically):

### Iopamidol (6 tubes):
- 3 tubes at 20mM: pH 6.2, 6.8, 7.4
- 3 tubes at 50mM: pH 6.2, 6.8, 7.4

### Creatine (6 tubes):
- 3 tubes at 20mM: pH 6.2, 6.8, 7.4
- 3 tubes at 50mM: pH 6.2, 6.8, 7.4

### Taurine (6 tubes):
- 3 tubes at 20mM: pH 6.2, 6.8, 7.4
- 3 tubes at 50mM: pH 6.2, 6.8, 7.4

### PLL (3 tubes):
- 0.1% w/v at pH 6.2, 6.8, 7.4

### PBS Blanks (3 tubes):
- PBS at pH 6.2, 6.4, 6.8

## 🚀 Quick Start Guide

### Option 1: B0 Correction Only (Faster, Recommended for Initial Analysis)

```matlab
% In MATLAB, navigate to the BM_sim_fit directory
cd /path/to/BM_sim_fit

% Run the B0-corrected analysis
phantom_CEST_analysis_B0corr
```

**You will be prompted to select:**
1. CEST Z-spectrum DICOM folder
2. B0 map DICOM folder (optional, but recommended)

**Processing time:** ~5-10 minutes

### Option 2: Full Correction (Most Accurate, Recommended for Final Results)

```matlab
% Run the full analysis with all corrections
phantom_CEST_analysis_FULL
```

**You will be prompted to select:**
1. CEST Z-spectrum DICOM folder
2. T1 map DICOM folder
3. T2 map DICOM folder
4. B1 map folder (NIfTI or DICOM)
5. B0 map DICOM folder

**Processing time:** ~15-30 minutes (depending on whether T1/T2 maps need calculation)

## 📁 Required Data Organization

Organize your DICOM data in separate folders:

```
Your_Data_Directory/
├── phantom_cest_v1/         # CEST Z-spectrum DICOMs
│   ├── IM_0001.dcm
│   ├── IM_0002.dcm
│   └── ...
├── phantom_t1/              # T1-weighted images (multiple TRs)
│   ├── IM_0001.dcm
│   └── ...
├── Phantom_t2/              # T2-weighted images (multiple TEs)
│   ├── IM_0001.dcm
│   └── ...
├── B1_map_dicom/            # B1 inhomogeneity map
│   ├── B1map.nii  OR
│   └── IM_0001.dcm
└── B0_map/                  # B0 field map
    ├── IM_0001.dcm
    └── ...
```

## 🔍 What the Scripts Do

### Automatic Processing Steps:

1. **Load CEST Data**
   - Reads all DICOM files
   - Extracts or prompts for PPM offsets
   - Creates Z-spectrum volume

2. **Load Correction Maps**
   - T1 relaxation map (if available)
   - T2 relaxation map (if available)
   - B1 inhomogeneity map (if available)
   - B0 field map (if available)

3. **Automatic Phantom Detection**
   - Detects phantom outline using Otsu thresholding
   - Identifies 24 tubes using connected component analysis
   - Orders tubes circularly around phantom center
   - **No manual ROI drawing required!**

4. **Apply Corrections**
   - B0: Shifts Z-spectra to correct frequency offsets
   - B1: Corrects for RF field inhomogeneity
   - T1: Accounts for longitudinal relaxation
   - T2: Accounts for transverse relaxation

5. **Calculate %CEST**
   - Computes MTR asymmetry for multiple offsets:
     - 1.9 ppm (Creatine)
     - 2.0 ppm (Creatine region)
     - 3.25 ppm (Taurine)
     - 3.5 ppm (Amine region)
     - 4.3 ppm (Iopamidol amide I)
     - 5.6 ppm (Iopamidol amide II)

6. **Estimate Exchange Rates (Kex)**
   - Uses simplified Bloch-McConnell model
   - Estimates Kex for each tube

7. **Generate Outputs**
   - Parametric maps (PNG images)
   - CSV file with all statistics
   - MATLAB .mat file with tube masks

## 📊 Output Files

### Version 1 (B0 correction):
```
CEST_results_B0corrected.csv           # Main results table
CEST_phantom_outline.png               # Phantom detection
CEST_tube_detection.png                # 24 tubes labeled
CEST_parametric_map_1.90ppm.png       # Creatine map
CEST_parametric_map_3.25ppm.png       # Taurine map
CEST_parametric_map_4.30ppm.png       # Iopamidol map
CEST_parametric_map_5.60ppm.png       # Iopamidol map
CEST_summary_heatmap.png               # All tubes × all offsets
CEST_tubeMasks.mat                     # Tube masks for reuse
```

### Version 2 (Full correction):
```
CEST_results_FULL_corrected.csv        # Comprehensive results
FULL_tube_detection.png                # Detection visualization
FULL_CEST_map_*.png                    # CEST maps for each offset
FULL_Kex_parametric_map.png            # Exchange rate map
FULL_tubeMasks.mat                     # Tube masks
T1map.mat                              # Calculated T1 map (cached)
T2map.mat                              # Calculated T2 map (cached)
```

## 📈 Understanding the CSV Output

The CSV file contains the following columns:

| Column | Description |
|--------|-------------|
| `tube_number` | Tube index (1-24) |
| `tube_label` | Chemical name, concentration, pH |
| `num_voxels` | Number of pixels in tube ROI |
| `T1_mean_ms` | Mean T1 relaxation time (ms) |
| `T2_mean_ms` | Mean T2 relaxation time (ms) |
| `B1_mean_percent` | Mean B1 field strength (%) |
| `B0_mean_ppm` | Mean B0 field offset (ppm) |
| `CEST_1.90ppm_mean` | %CEST at 1.9 ppm (Creatine) |
| `CEST_1.90ppm_std` | Standard deviation |
| `CEST_3.25ppm_mean` | %CEST at 3.25 ppm (Taurine) |
| `CEST_4.30ppm_mean` | %CEST at 4.3 ppm (Iopamidol) |
| ... | (continues for all offsets) |
| `Kex_Hz` | Estimated exchange rate (Hz) |

## ⚙️ Customization Options

### Modify Target CEST Offsets

Edit the `target_offsets` array in the scripts:

```matlab
% Default offsets
target_offsets = [1.9, 2.0, 3.25, 3.5, 4.3, 5.6];

% To add or change offsets:
target_offsets = [1.9, 2.0, 3.25, 3.5, 4.3, 5.6, 2.75];  % Add 2.75 ppm
```

### Modify Tube Labels

If your tube order is different, edit the `tube_labels` cell array at the top of the script.

### Adjust Colormap Ranges

Change the colormap limits for parametric maps:

```matlab
% Default CEST range
caxis([0 10]);  % 0-10% CEST

% To change:
caxis([0 15]);  % For higher CEST effects
```

### Change Field Strength

If you're not using a 3T scanner, modify:

```matlab
B0_field = 3.0;  % Tesla

% For 7T:
B0_field = 7.0;
```

## 🔧 Troubleshooting

### Problem: "Expected 24 tubes but detected X"

**Solution:**
- Check image quality and contrast
- Adjust Gaussian smoothing parameter:
  ```matlab
  bgd = imgaussfilt(S0_image, 4);  % Try values 3-6
  ```
- Adjust minimum tube size threshold:
  ```matlab
  indexOfSmall = find(numOfPixels < 10);  % Try 5-20
  ```

### Problem: "PPM offsets not found in DICOM headers"

**Solution:**
- Choose "Load from file" and provide a text/CSV file with offsets
- Or choose "Manual Entry" and input offsets for each file
- Or use "Default Range" for testing (-5 to +5 ppm)

### Problem: "No B0/B1/T1/T2 map found"

**Solution:**
- The script will continue without that correction
- You can skip optional maps and still get results
- For publication-quality results, acquire all correction maps

### Problem: Slow processing with high-resolution images

**Solution:**
- Scripts automatically use voxel-wise processing
- Consider downsampling if needed:
  ```matlab
  zspecVolume = imresize3(zspecVolume, 0.5);  % 50% resolution
  ```

### Problem: Tube labels don't match physical layout

**Solution:**
- The automatic detection orders tubes clockwise from an arbitrary starting point
- After running the script once, note which tube is which from the output images
- Manually reorder the `tube_labels` array to match

## 📝 Tips for Best Results

1. **Use high SNR reference image**: The S0 image (near 0 ppm) should have good contrast

2. **Acquire B0 map**: B0 correction is the most important for accurate CEST

3. **Consistent slice positioning**: Ensure all maps (T1/T2/B1/B0/CEST) are from the same slice

4. **Check tube detection**: Always review `CEST_tube_detection.png` to verify correct detection

5. **Reuse tube masks**: Once detected correctly, save and reuse `tubeMasks.mat`:
   ```matlab
   load('CEST_tubeMasks.mat');  % Load pre-computed masks
   ```

6. **Validate results**: Compare CEST values to expected ranges:
   - Creatine @ 1.9 ppm: 2-5% CEST
   - Iopamidol @ 4.3 ppm: 5-15% CEST (concentration dependent)
   - PBS blanks: <1% CEST

## 📚 Expected CEST Values (Reference)

| Chemical | Concentration | pH | Expected CEST @ 3T | Exchange Regime |
|----------|--------------|-----|-------------------|-----------------|
| Iopamidol | 20mM | 6.2 | 5-8% @ 4.3ppm | Slow |
| Iopamidol | 50mM | 6.2 | 10-15% @ 4.3ppm | Slow |
| Creatine | 20mM | 6.8 | 2-4% @ 1.9ppm | Intermediate |
| Creatine | 50mM | 6.8 | 4-7% @ 1.9ppm | Intermediate |
| Taurine | 20mM | 6.8 | 1-3% @ 3.25ppm | Fast |
| Taurine | 50mM | 6.8 | 3-5% @ 3.25ppm | Fast |
| PLL | 0.1% | 6.8 | 3-6% @ 3.6ppm | Intermediate |
| PBS | - | - | <1% | - |

**Note:** Values vary with saturation power, duration, and B0 field strength.

## 🔬 Advanced: Bloch-McConnell Fitting

For more accurate Kex estimation, you can use the full BM fitting tools in the `optimisation/` directory:

```matlab
% After running the analysis script:
load('CEST_tubeMasks.mat');

% Extract Z-spectrum for tube 1
tube1_mask = tubeMasks(:,:,1);
tube1_spectrum = mean(zspecVolume(repmat(tube1_mask, [1 1 numFiles])));

% Prepare for BM fitting
Zlab.x = ppmOffsets;
Zlab.y = tube1_spectrum / max(tube1_spectrum);  % Normalize

% Run multiZfit (requires proper Sim structure setup)
% See doc/BM_Documentation.docx for details
```

## 📞 Support & References

### Key Functions Used:
- `t1fitting_VTR.m` - Variable TR T1 fitting
- `t2fitting.m` - Multi-echo T2 fitting
- `calcMTRmap.m` - MTR asymmetry calculation
- `multiZfit.m` - Bloch-McConnell parameter fitting

### Useful MATLAB Commands:
```matlab
% View tube detection results
load('CEST_tubeMasks.mat');
figure; imagesc(sum(tubeMasks, 3)); colorbar;

% Load and inspect CSV results
T = readtable('CEST_results_FULL_corrected.csv');
summary(T)

% Plot Z-spectrum for specific tube
tube_idx = 5;
plot(ppmOffsets, squeeze(mean(mean(zspecVolume .* tubeMasks(:,:,tube_idx), 1), 2)));
xlabel('Frequency offset (ppm)');
ylabel('Normalized signal');
title(sprintf('Z-spectrum: %s', tube_labels{tube_idx}));
```

## 🎯 Workflow Recommendation

For your first analysis:

1. **Run T1 and T2 mapping scripts separately first**:
   ```matlab
   phantom_t1      % Generates T1map.mat
   phantom_T2      % Generates T2map.mat
   ```

2. **Run B0-corrected version for quick check**:
   ```matlab
   phantom_CEST_analysis_B0corr
   ```

3. **Review outputs and verify tube detection**

4. **Run full correction version for final analysis**:
   ```matlab
   phantom_CEST_analysis_FULL
   ```

5. **Export final CSV and parametric maps**

---

## 📄 Citation

If you use this pipeline in your research, please cite:

- The original BM_sim_fit repository
- Any relevant CEST methodology papers

## ⚖️ License

Follows the license of the parent BM_sim_fit repository (see LICENSE.md).

---

**Last Updated:** 2025-11-05
**Version:** 1.0
**Contact:** See repository issues for support
