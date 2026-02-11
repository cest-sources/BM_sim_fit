# Bloch-McConnell CEST Analysis Pipeline

## New Files for BMC-Based Analysis

### 1. `phantom_CEST_BMC_analysis.m` - Main Analysis Script

**What it does:**
- **Manual ROI drawing** for 24 tubes (in correct order)
- WASSR-based B0 field correction
- Bloch-McConnell-based analysis (replaces Lorentzian fitting)
- Kex (exchange rate) extraction
- Concentration estimation
- Parametric map generation

**Key Features:**
- ✅ No Lorentzian fitting
- ✅ BMC simulation for %CEST calculation
- ✅ Kex extraction from CEST effects
- ✅ Concentration estimation per tube
- ✅ WASSR B0 correction
- ✅ **Manual ROI drawing** (draw tubes 1-24 in exact order)
- ✅ Interactive drawing with real-time preview
- ✅ Auto-saves backups every 6 tubes

### 2. `reorder_tubes_to_match_layout.m` - Tube Reordering Helper

**When to use:**
- **Not needed** with manual ROI drawing (you draw tubes in correct order)
- Only useful if you used automatic detection in a previous version
- Can reorder pre-existing tube masks to match your phantom layout

## Quick Start Guide

### Step 1: Run Main Analysis

```matlab
cd /path/to/BM_sim_fit
phantom_CEST_BMC_analysis
```

### Step 2: Select ROI Drawing Method

When prompted, choose your preferred drawing method:
- **Option 1 (Recommended)**: Circle ROI - Best for round tubes
- **Option 2**: Polygon ROI - For irregular shapes
- **Option 3**: Freehand ROI - For complex shapes

### Step 3: Draw ROIs for All 24 Tubes

**Have your shots.png image open for reference!**

The script will guide you through drawing ROIs for each tube in order:

1. **Tube 1**: Iopamidol 20mM pH 6.2
2. **Tube 2**: Iopamidol 20mM pH 6.8
3. **Tube 3**: Iopamidol 20mM pH 7.4
...and so on through Tube 24

**For each tube:**
- The screen shows the tube number and chemical name
- Previously drawn tubes are displayed with numbers
- Draw the ROI around the current tube
- Double-click inside to finish
- Confirm or redraw as needed
- Progress auto-saves every 6 tubes

### Step 4: Verify and Continue

After drawing all ROIs:
- Script shows overview with all tubes numbered
- Verify numbering matches your shots.png
- Analysis automatically continues with BMC fitting

### Step 5: Review Results

Check the generated files:
- `BMC_CEST_results.csv` - Kex, concentration, %CEST values
- `BMC_parametric_maps.png` - All maps
- `MANUAL_tube_numbering.png` - Your ROIs with numbers

## Your Phantom Layout

Based on your shots.png image:

| Tube Numbers | Chemical | Concentration | pH Values |
|--------------|----------|---------------|-----------|
| 1-3 | Iopamidol | 20 mM | 6.2, 6.8, 7.4 |
| 4-6 | Iopamidol | 50 mM | 6.2, 6.8, 7.4 |
| 7-9 | Creatine | 20 mM | 6.2, 6.8, 7.4 |
| 10-12 | Creatine | 50 mM | 6.2, 6.8, 7.4 |
| 13-15 | Taurine | 20 mM | 6.2, 6.8, 7.4 |
| 16-18 | Taurine | 50 mM | 6.2, 6.8, 7.4 |
| 19-21 | PLL | 0.1% w/v | 6.2, 6.8, 7.4 |
| 22-24 | PBS (blank) | - | 6.2, 6.8, 7.4 |

## Output Files

After running the analysis, you'll get:

1. **BMC_CEST_results.csv** - Main results table with:
   - Tube numbers and labels
   - Kex values (Hz)
   - Concentration estimates (mM)
   - %CEST at multiple offsets (1.9, 3.5, 4.3 ppm)
   - B0 field values per tube

2. **BMC_parametric_maps.png** - Comprehensive visualization:
   - Kex map
   - %CEST maps at 4.3 ppm (Iopamidol)
   - %CEST maps at 1.9 ppm (Creatine)
   - %CEST maps at 3.5 ppm (Amide)
   - B0 field map
   - Kex bar chart by tube

3. **DETECTED_tube_numbering.png** - Tube detection verification image

4. **BMC_tubeMasks.mat** - Tube masks for reuse

5. **BMC_CEST_workspace.mat** - Full MATLAB workspace

## CSV Output Format

| Column | Description |
|--------|-------------|
| `tube_number` | 1-24 |
| `tube_label` | Chemical, concentration, pH |
| `num_voxels` | Number of pixels in tube ROI |
| `Kex_Hz` | Exchange rate (Hz) |
| `concentration_mM` | Estimated concentration (mM) |
| `chemical_shift_ppm` | Expected chemical shift |
| `CEST_at_3_5ppm` | %CEST at 3.5 ppm (amide) |
| `CEST_at_4_3ppm` | %CEST at 4.3 ppm (Iopamidol) |
| `CEST_at_1_9ppm` | %CEST at 1.9 ppm (Creatine) |
| `B0_mean_ppm` | Mean B0 field offset |

## Differences from Previous Scripts

### ❌ Removed:
- Lorentzian multi-pool fitting
- Manual ROI selection for single tube
- Generic offset handling

### ✅ Added:
- Bloch-McConnell-based Kex extraction
- Concentration estimation
- Automatic processing of all 24 tubes
- Tube numbering verification system
- Chemical-specific CEST measurements
- Predefined offsets for your acquisition protocol

## Kex Calculation Method

The script uses a simplified Bloch-McConnell approach:

```
Kex ≈ 2π × Δω × (MTRasym% / f_b)
```

Where:
- Δω = chemical shift frequency (rad/s)
- MTRasym = measured CEST effect (%)
- f_b = bound pool fraction (~0.5% assumed)

For more accurate Kex values, you can integrate with `multiZfit.m` from the optimisation directory.

## Troubleshooting

### Issue: "Accidentally drew wrong tube"

**Solution:**
- Click "Redraw" when prompted
- Or click "Skip" and redraw it later manually

### Issue: "Forgot which tube I'm on"

**Solution:**
- Check the title bar - shows "TUBE X/24: Chemical name pH"
- Previously drawn tubes are shown with numbers
- If unsure, click "Skip" and fix later

### Issue: "Tube numbering doesn't match shots.png after drawing"

**Solution:**
- You likely drew tubes in wrong order
- Delete `BMC_tubeMasks.mat` and rerun script
- OR use `reorder_tubes_to_match_layout.m` to fix the order

### Issue: "Kex values seem too high/low"

**Solution:**
- Check B0 correction quality (view B0_map)
- Verify S0 normalization
- Adjust assumed bound pool fraction (currently 0.5%)

### Issue: "Want to use full multiZfit instead of simplified model"

**Solution:**
- Uncomment multiZfit integration sections in the BMC script
- Ensure all BM simulation functions are in MATLAB path
- May need to adjust Sim structure parameters

## Advanced Usage

### Integrate Full Bloch-McConnell Fitting

To use the complete `multiZfit` function for more accurate Kex:

1. Ensure optimisation folder is in path
2. Set up Sim structure with appropriate parameters:

```matlab
Sim = init_Sim();
Sim.B0 = 3.0;  % Tesla
Sim.FREQ = [0; expected_shift_ppm];  % Water + CEST pool
Sim.T1 = [2.0; 1.5];  % T1 values (s)
Sim.T2 = [0.08; 0.01];  % T2 values (s)
% ... continue setup ...

FIT = multiZfit(tube_zspec, ppmOffsets, Sim);
Kex_fitted = FIT.results.k_BA;
```

3. Replace simplified Kex calculation in Step 7

## Support

For issues or questions about:
- **Tube detection**: Check Step 5 in main script
- **B0 correction**: Check Step 3 and 6
- **BMC fitting**: Check Step 7
- **Tube reordering**: Use helper script

## Version History

- **v1.0** (2025-11-05): Initial BMC-based analysis pipeline
  - Removed Lorentzian fitting
  - Added BMC Kex extraction
  - Added tube numbering verification
  - Integrated WASSR B0 correction

---

**Last Updated:** 2025-11-05
**Requires:** MATLAB R2019b or later with Image Processing Toolbox
