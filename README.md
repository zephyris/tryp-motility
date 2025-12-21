
# tryp-motility
ImageJ macros for analysis of trypanosomatid (_Trypanosoma brucei_, _Leishmania spp._, etc.) from the paper: Wheeler RJ (2017) “Use of chiral cell shape to ensure highly directional swimming in trypanosomes” PLoS Comput. Biol. 13(1):e1005353 [doi:10.1371/journal.pcbi.1005353](https://doi.org/10.1371/journal.pcbi.1005353).

## Requirements and installation
This requires any recent version of [ImageJ](https://imagej.net/ij/) or [Fiji](https://fiji.sc/), and should work on any ImageJ-supported OS. These scripts were tested with ImageJ 1.53 on Windows 11.

Install using `Plugins > Macros > Install` and selecting `ParticleTrackingv4.ijm`. This adds menu options to the `Plugins` menu and should load essentially instantly.

## Usage
For basic usage, open a multi-slice TIFF video of swimming cells. See [method readme](method.md) for suggested sample preparation and video capture details.

Configure detection and tracking using:
* `Particle Detection Settings [q]`, with a detection preview using `Preview Particle Detection [w]`.
* `Particle Track Analysis Settings [e]`

Carry out an analysis using:
* `Particle Motility Analysis [a]`, and field mean drift can be subtracted using `Subtract Drift From Particle Tracks`.
Runtime heavily depends on the size of the input video. For the example data (1024x768x64) it takes ~15 seconds.

Results are output in the log window. This is a table where each row represents a cell track, with self-explanatory column headings. Additional visualisations can be plotted using:
* `Particle Track Plotter [z]` (the motion of all particles in the field of view)
* `Display Track [x]` (make a line selection for the selected cell number)
* `Extract Cell [c]` (extract a cropped view of the selected cell number)

## Batch usage
All files in a directory can be analysed using:
* `Batch Analyse Folder`

This expects to find `.tif` or `.czi` (Zeiss proprietory) images. Loading `.czi` images requires the [Bio-Formats plugin](https://imagej.net/formats/bio-formats).

Output data is a text file containing the data table and images of the particle tracks, with and without drift correction.

`QuickSummariser.ijm` carries out a track length-weighted average of data in all standard output tables in a directory.

## Very batch usage
For parallel headless analysis (Windows only), copy `ParticleTrackingv4.bat` and `ParticleTrackingv4.ijm` to the root directory of an ImageJ installation, and copy `ParticleTrackingBatchv4.bat` to a data directory. Ensure all paths are configured in these scripts, then run `ParticleTrackingBatchv4.bat` for parallel headless analysis of all files in the directory.

## Example data
A short cropped example video (`LmxWT_5fps_crop-short.tif`) is provided, along with the standard outputs from batch analysis in `ExampleData`.
This video is at 5 frames per second at a magnification of 1.5384 pixels per micron.
