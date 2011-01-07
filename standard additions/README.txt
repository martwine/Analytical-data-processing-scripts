data files

Desc: text string describing data. Single word description.
rowflag: 1 - data included in calculation by default; 0 - data not included by default
set:	denotes set of associated measurements i.e. one standard addition calibration, baseline measurements, uncalibrated unknowns
type: 	positive integer - volume of stock standard added in uL
	-999 unknown requiring calibration by preceeding std addition calibration in file
	-555 baseline
	-666 alternative baseline (e.g. for different approaches to background fluorescence in NH4 analysis)
	-333 Tracer 

1,2,3,4,5,6,7,8,9,10 : values of repeat measurements (only 1 required)

method files are a single row of data containing the folling columns

stock_concn	
sample_volume	
baseline_method 
tracer_name	
tracer_type
metadata


