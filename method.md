This is a recommended protocol for generating videos for analysis using this code.

## Method
Swimming can vary with culture density and health. Make sure that the cultures are growing nicely - log phase growth with daily subculture to a defined density. Make sure any controls are grown in parallel to experimental samples.
Make sure timing is correct that cultures are sampled at a good density (to measure many cells) in log phase. For example, for T. brucei,  aim for 5x10<sup>6</sup> cells/ml (eg. end of day subculture to 1x10>sup>6</sup>, then analyse the following morning).

Samples are best analysed in a deep (compared to focal depth) chamber to avoid surface effects.
Adhesive 1 x 1 cm, 0.25 mm deep "Gene Frames" by Thermo (product code 11560294) work well.
Prepare a standard glass slide by attaching the Gene Frame sticker to make the well/chamber to hold the cells.

Take 5 μl culture medium and put in the middle of the square well/chamber, and immediately adhere coverslip on top (to avoid evaporation).
To avoid liquid flow due to capillary action, aim to keep the cells as a droplet between the slide and coverslip not touching the sides of the chamber/well. Imagine a circle drop in the middle of a square chamber.
If the droplet touches the side of the chamber then capillary action will tend to make the liquid flow rapidly compared to the swimming speed of the cells.

Transfer the sample to the microscope. Some temperature shock is unavoidable (the sample is so small that it immediately cools to room temperature), so consistent timing of handling ensures consistency.
Focus in the middle (x, y and z) of the droplet in the chamber. It is best, if possible, to use the microscope stage height measurement to guide focusing in the vertical centre of the sample.
Capture the video.

There is one key sample handling detail to consider. Immotile cells will gradually settle onto the slide over time, so might not be represented in the video.
The best way to compensate for this is rigorous sample handling. Eg. Prepare the slide with coverslip, then turn upside down (so any immotile cells are gradually settling onto the coverslip). Carry to the microscope upside down, then flip to put on the microscope stage. Wait a standard length of time (eg. 30s), then take the video.

Keep the samples growing in culture for three or four days with regular subculture, and take a video each day, to capture good biological replicates.

## Microscopy
Recoomended settings are:
* 512 frame video
* 5 frames per second
* 10x magnification
* Darkfield microscopy
Make sure that the frame rate is a true frame rate, not an interval, ie. 100ms exposure time and 400ms wait is correct, not 100ms exposure time and 500ms wait.
Darkfield microscopy often possible on microscopes with no special configuration, using a condenser Ph3 ring and standard 10x objective.

An uncompressed (eg. TIFF) video is much better than a compressed (eg. AVI, MP4, MPG, etc.) video. Video compression is not well-optimised for small moving particles and adds artefacts.

## Output
Analysis outputs a data file with one row per cell track. When considering this data, it is best to weight by track length in seconds to avoid biasing towards swimming behaviour which tend to give short tracks.

The simplest analysis is mean swimming speed from each of three or four videos. These those three or four data points can be compared using t-test (as central limit theorem is met) to a control sample.
