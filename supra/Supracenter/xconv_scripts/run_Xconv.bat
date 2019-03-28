@echo off
rem Batch file to convert Met Office .pp files to NetCDF format
rem Karen Amoudry, 18 June 2013, National Oceanography Centre, Liverpool
setlocal EnableDelayedExpansion
call :treeProcess

:treeProcess
	for %%f in (*.pp) do (
		set outname=%%~nf.nc
		set "cmd_str=convsh1.94.exe < conv2nc.tcl -i %%~ff -o !outname!"
		echo !cmd_str!
		)


echo ^@echo off
for /D %%d in (*) do (
	cd %%d
	call :treeProcess
	cd ..
)
exit /b