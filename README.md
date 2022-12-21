# Analysis pipeline for "Test-Retest Reliability of the Human Connectome: An OPM-MEG study" (In submission 2022)

Lukas Rier[1], Sebastian Michelmann[2], Harrison Ritz[2], Vishal Shah[3], Ryan M. Hill[1,4], James Osborne[3], Cody Doyle[3], Niall Holmes[1,4], Richard Bowtell[1], Matthew J. Brookes[1,4], Kenneth A. Norman[2], Uri Hasson[2], Jonathan D. Cohen[2] and Elena Boto[1,4]
 
[1] Sir Peter Mansfield Imaging Centre, School of Physics and Astronomy, University of Nottingham, University Park, Nottingham, NG7 2RD, UK.

[2] Princeton Neuroscience Institute, Washington Road, Princeton, NJ 08544, USA.

[3] QuSpin Inc. 331 South 104th Street, Suite 130, Louisville, Colorado, 80027, USA.

[4] Cerca Magnetics Limited, 2 Castlebridge Office Village, Kirtley Drive, Nottingham, NG7 1LD, Nottingham, UK.

-------

All data used can be found at https://doi.org/10.5072/zenodo.1134455.
Create a new directory - your project_directory - containing subdirectories 'scripts' and 'data'

Download FieldTrip version 20199212 at https://www.fieldtriptoolbox.org/download.php and save the toolbox to '/project_directory/scripts'.

Download and extract the zip folder to '/project_directory/data' and clone this repository to '/project_directory/scripts' to get the following file structure:
   
|project_directory   
|---|--> data   
|------|--> derivatives   
|------|--> sub-001   
|------|--> sub-002   
|------|--> ...   
|------|--> sub-010   
|------participants.tsv   
|---|--> scripts   
|------|--> Beamformer   
|------|--> BrainPlots   
|------|--> fieldtrip-20199212   
|------|--> gifti-1.8   
|------...   

Run_VEs_AEC.m will generate all necessary outputs to produce results presented in "Test-Retest Reliability of the Human Connectome: An OPM-MEG study".
Figure_Spectral_Power.m, Figure_individual_connectomes.m, Figure_Conn_fingerprinting.m,Figure_average_connectivity.m reproduce all figure elements
Edit the scripts above to include the path to your project directory in ```project_dir``` before running the analyses.
