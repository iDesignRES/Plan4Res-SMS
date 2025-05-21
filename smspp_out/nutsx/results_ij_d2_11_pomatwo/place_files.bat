cls
@pushd %~dp0

move /Y ActivePowerOUT.csv ActivePower/.
move /Y DemandOUT.csv Demand/.
for %%f in (Flows*.csv) do (
    echo "name %%f"
    move /Y %%f Flows/.
)
move /Y MaxPowerOUT.csv MaxPower/.
move /Y PrimaryOUT.csv Primary/.
move /Y SecondaryOUT.csv Secondary/.
move /Y VolumeOut.csv Volume/.
for %%f in (MarginalCost*.csv) do (
    echo "name %%f"
    move /Y %%f MarginalCosts/.
)
@popd
rem pause


