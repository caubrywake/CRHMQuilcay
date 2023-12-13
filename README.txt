The file structure for the CRHM Quilcayhuanca project is as follow:

Main Folder: CRHMcuchi

Subfolder and files:

- CRHM: version 
	* Output (Created from R script, running CRHM simulations)
		Noice_RunFlow.txt
		NoLeakage_RunFlow.txt
		RunET.txt
		RunFlow.txt
		RunForcings.txt
		RunRouting.txt 
		NoGW_RunFlow.csv (Output from CRHM processed to datetime table)
		NoGW_RunFlow_exportedfromCRHM.obs (this was exported manually from CRHM, so it has a slighlty different format)
		* scenario
			(90 .txt files from CRHM output)
	* prjfile
		quilcay_20230823.prj
		quilcay_20230823_noice.prj
		quilcay_20230823_noleakage.prj
		quilcay_20230823_nogw.prj
		quilcay_scenario_20230823.prj
		setvariables_Basinflow.prj
		setvariables_Forcings.prj
		setvariables_Scenarios.prj
		setvariables_ET.prj
		setvariables_Routing.prj
- scripts
	AssembleObs_supplementary.m
	CalculatingLascarTempGradient_Supplementary.m
	RunningCRHM_Cuchi_20230823.R
	ImportingCRHMoutput.m
	PlottingWells.m
	ModelEval_streamflow.m
	Scenario.m
	LikelyScenario_timeseries.m

- fig (folder where figures are saved)
 	modeleval (empty folder to receive figures from code)
	scenario (empty folder to receive figures from code)
	obs (empty folder to receive figures from code)
	hru (with maps of the HRU)

- data
	*Q
		cda_lev-q-t.csv (streamflow data)
	*well	
		(all .csv)
	*processed
		CuchiObs_2014_2020.obs (CRHM obs file)
		
- functions
	nashsutcliffe.m
	klinggupta.m