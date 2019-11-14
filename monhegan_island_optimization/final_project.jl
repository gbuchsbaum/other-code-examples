# Model for Gabriel Buchsbaum's project on microgrid design for Monhegan Island

####################################
######### Initialize tools #########
####################################
using JuMP
using AmplNLWriter
using Cbc
using DataFrames
using CSV
using Statistics

##########################################
####### Read and process load data #######
##########################################

# Currently using historical data from ISO New England, rather than
# location-specific data. Data is in MWh
LoadData = CSV.read("ISO_data.csv", header=5, datarow=7)

# Identify number of timesteps (i.e. number of hours in a year)
nTimes = length(LoadData[1])

# Average summer load of the island [kW]
SummerAverage = 80
SummerFactor = SummerAverage/Statistics.mean(LoadData[4][5233:5257])
WinterAverage = 27
WinterFactor = WinterAverage/Statistics.mean(LoadData[4][8425:8449])
AverageFactor = 310000 / sum(LoadData[4])
# Convert the ISO data to the correct scale [kW]
Load = zeros(nTimes)
for t=1:nTimes
    if t >= 3337 && t <= 6049
        Load[t] = LoadData[4][t] * SummerFactor
    else
        Load[t] = LoadData[4][t] * WinterFactor
    end
end
# Planning Reserve Margin (amount capacity needs to exceed peak load for reliability)
PRM = 0.15

###########################################
####### Read and process solar data #######
###########################################

SolarFile = "1397487_43.77_-69.3_2013.csv"

# Read solar data file. File is dowloaded using NREL's System Advisor Model
SolarData = CSV.read(SolarFile, header=3)

# Convert days of the month to days of the year by increasing the day number
# every time a new day begins
Day = zeros(Int64,nTimes)
Day[1] = 1
for i = 2:nTimes
    if SolarData[3][i] == SolarData[3][i-1]
        Day[i] = Day[i-1]
    else
        Day[i] = Day[i-1]+1
    end
end

# Make an array of the clock hour used for solar data
Hour = zeros(nTimes)
for i = 1:nTimes
    Hour[i] = SolarData[4][i] + SolarData[5][i]/60
end

# Identify the three main measures of irradiance from the file.
# Sunlight traveling directly to a perpendicular surface [W/m2]
DirectNormal = SolarData[7]
# Sunlight scattered by the atmosphere, coming equally from all directions [W/m2]
DiffuseHorizontal = SolarData[6]
# Total sunlight reaching a horizontal surface [W/m2]
GlobalHorizontal = SolarData[8]

# Fraction of light reflected
Albedo = SolarData[10]

# Read the header of the file to identify the latitude and longitude
SolarData2 = CSV.read(SolarFile; rows=5)
# Due to future calclulation needs, latitude is stored in radians,
# while longitude is stored in degrees.
Latitude = deg2rad(parse(Float64,SolarData2[6][1]))
Longitude = parse(Float64,SolarData2[7][1])

# Identify location of sun at all times of the year.
# Since Julia uses radians for trigonometric functions, all angles are converted
# to radians
AltitudeAngle = zeros(nTimes)
AzimuthAngle = zeros(nTimes)
for i = 1:nTimes
    # Determine the offset between clock time and solar time, based on the time
    # of year and location in the time zone
    B = 2 * pi / 364 * (Day[i] - 81)
    EquationOfTime = 9.87 * sin(2*B) - 7.53*sin(B) - 1.5 * sin(B)
    SolarTime = Hour[i] + 4/60 * (Longitude + 75) + EquationOfTime/60
    # Find the angle based on the time until solar noon
    HourAngle = deg2rad(15 * (12 - SolarTime))
    # Find the Earth's tilt relative to the sun
    Declination = deg2rad(23.45 * sin(2*pi*(Day[i]-81)/365))
    # Find the height of the sun above the horizon
    AltitudeAngle[i] = asin(cos(Latitude) * cos(Declination) * cos(HourAngle) + sin(Latitude) * sin(Declination))
    # Find the east-west position of the sun
    if cos(HourAngle) >= tan(Declination) / tan(Latitude)
        AzimuthAngle[i] = asin(cos(Declination) * sin(HourAngle) / cos(AltitudeAngle[i]))
    else
        AzimuthAngle[i] = pi - asin(cos(Declination) * sin(HourAngle) / cos(AltitudeAngle[i]))
    end
end

PATH_TO_SOLVERS = ENV["ERE291_SOLVERS"]
ampl_solver = "snopt"
m1 = Model(solver=AmplNLSolver(joinpath(PATH_TO_SOLVERS,ampl_solver),["outlev=2"]))

# The solar collector tilt
@variable(m1, CollectorTilt >= 0,start=deg2rad(20))

# Direction collectors are facing
@variable(m1, CollectorAzimuth, start= 0)

# Identify the portion of the rated capacity that is available at each time

# Find the angle between the collector and sun
@NLexpression(m1, cosTheta[t=1:nTimes], cos(AltitudeAngle[t]) * cos(AzimuthAngle[t] - CollectorAzimuth) * sin(CollectorTilt) + sin(AltitudeAngle[t]) * cos(CollectorTilt))
# Find the direct normal radiation in the plane of the panels [W/m2]
@NLexpression(m1, DirectNormalGeneration[t=1:nTimes], DirectNormal[t] * cosTheta[t])
# Find the diffuse radiation seen by the panels (i.e. portion of sky seen) [W/m2]
@NLexpression(m1, DiffuseHorizontalGeneration[t=1:nTimes], DiffuseHorizontal[t] * (1 + cos(CollectorTilt))/2)
# Find the reflected radiation seen by the panels [W/m2]
@NLexpression(m1, ReflectedGeneration[t=1:nTimes], Albedo[t]*GlobalHorizontal[t] * (1 - cos(CollectorTilt))/2)
# Find the portion of rated capacity generated at that time, assuming that
# rated capacity is at 1 sun or 1000 W/m2 of sunlight
@NLexpression(m1, SolarGeneration[t=1:nTimes], (DirectNormalGeneration[t] + DiffuseHorizontalGeneration[t]+ReflectedGeneration[t])/1000)

# Optimize sum of products of generation and load
@NLobjective(m1, Max, sum(SolarGeneration[t]*Load[t] for t=1:nTimes))

# angles must be within reasonable ranges
@constraint(m1, CollectorTilt <= pi/2)
@constraint(m1, -pi/2 <= CollectorAzimuth <= pi/2)

solve(m1)

println("PV Tilt Angle: ",rad2deg(getvalue(CollectorTilt)), "°")
println("PV Azimuth Angle: ", rad2deg(getvalue(CollectorAzimuth)), "°")

# Array parameters. Consider adding these as variables to optimize,
# although they would require nonlinear optimization
# Tilt of collectors
CollectorTilt = getvalue(CollectorTilt)
# Direction collectors are facing
CollectorAzimuth = getvalue(CollectorAzimuth)
# Portion of power lost in conversion from solar cell to AC power
SolarLosses = 0.14

# Solar costs. These numbers are from Lazard, based on fixed community-scale solar
# Upfront capital cost for installation [$/kW]
FixedCapitalCost = 2500.00
# Annual cost for maintenance and other needs [$/kW/year]
FixedAnnualCost = 16.00
# Price increase from switching to single-axis tracking
TrackingCharge = 1.32
# Interest rate
LoanInterest = 0.05
# Lifetime of solar project [years]
SolarLife = 25
# Capital recovery factor
CRFsolar = LoanInterest * (1+LoanInterest)^SolarLife / ((1+LoanInterest)^SolarLife - 1)

# Identify the portion of the rated capacity that is available at each time
FixedGeneration = zeros(nTimes)
for i = 1:nTimes
    # Find the angle between the collector and sun
    cosTheta = cos(AltitudeAngle[i]) * cos(AzimuthAngle[i] - CollectorAzimuth) * sin(CollectorTilt) + sin(AltitudeAngle[i]) * cos(CollectorTilt)
    # Find the direct normal radiation in the plane of the panels [W/m2]
    DirectNormalGeneration = DirectNormal[i] * cosTheta
    # Find the diffuse radiation seen by the panels (i.e. portion of sky seen) [W/m2]
    DiffuseHorizontalGeneration = DiffuseHorizontal[i] * (1 + cos(CollectorTilt))/2
    # Find the reflected radiation seen by the panels [W/m2]
    ReflectedGeneration = Albedo[i]*GlobalHorizontal[i] * (1 - cos(CollectorTilt))/2
    # Find the portion of rated capacity generated at that time, assuming that
    # rated capacity is at 1 sun or 1000 W/m2 of sunlight
    FixedGeneration[i] = (1-SolarLosses)*(DirectNormalGeneration + DiffuseHorizontalGeneration+ReflectedGeneration)/1000
end

# Identify the portion of rated capacity for single-axis tracking at each time
TrackingGeneration = zeros(nTimes)
for i = 1:nTimes
    # Find the angle between the collector and the sun
    cosTheta = sqrt(1-(cos(AltitudeAngle[i]) * cos(AzimuthAngle[i]))^2)
    # Find the direct normal radiation in the plane of the panels [W/m2]
    DirectNormalGeneration = DirectNormal[i] * cosTheta
    # Find the diffuse radiation seen by the panels (i.e. portion of sky seen) [W/m2]
    DiffuseHorizontalGeneration = DiffuseHorizontal[i] * (1 + sin(AltitudeAngle[i]/cosTheta))/2
    # Find the reflected radiation seen by the panels [W/m2]
    ReflectedGeneration = Albedo[i]*GlobalHorizontal[i] * (1 + sin(AltitudeAngle[i]/cosTheta))/2
    # Find the portion of rated capacity generated at that time, assuming that
    # rated capacity is at 1 sun or 1000 W/m2 of sunlight
    TrackingGeneration[i] = (1-SolarLosses)*(DirectNormalGeneration + DiffuseHorizontalGeneration+ReflectedGeneration)/1000
end

##################################
####### Battery parameters #######
##################################

# The following are from (Lu, 2009) (Makarov 2008) (EPRI 2003)
# The cost of battery capacity [$/kWh]
CostCapacity = 250
# The cycle lifetime of the battery [kWh lifetime stored into battery per kWh of capacity]
CycleLife = 3000
# The round-trip losses [kWh/kWh stored]
RoundTripLosses = 0.0645

# The battery lifetime [years]
BatteryLifetime = 20
CRFbattery = CRFwind = LoanInterest * (1+LoanInterest)^BatteryLifetime / ((1+LoanInterest)^BatteryLifetime - 1)

##########################################
####### Read and process wind data #######
##########################################

# Read wind data file at 80m. Data is from NREL's SAM
WindData = CSV.read("lat43.76_lon-69.32__2013_40m.csv", header=3,datarow=6)
# Extract wind speed data [m/s]
Speed=WindData[3]

# Array to hold the name of the files storing wind turbine data
TurbineFiles = ["Wind_Turbine_Industries_Corp_Jacobs_31-20.csv","Pitchwind_30.csv","Ergycon_Ely50.csv","Deimos_60.csv","Northern_Power_Northwind_100.csv"]

# Descriptions of turbine options
TurbineOpts = ["20 kW","30 kW","50 kW","60 kW","100 kW"]

# Number of turbines being considered
nTurbineOpts = length(TurbineOpts)

# Array of DataFrames to hold hourly results
ResultsDF = Array{DataFrame}(undef, 2, nTurbineOpts)

# Array to hold total cost results
CostResults = zeros(2, nTurbineOpts)

# Array to hold amount of solar panels in solutions
kWSolarResults = zeros(2, nTurbineOpts)

# Array to hold number of wind turbines in solutions
nTurbinesResults = zeros(2, nTurbineOpts)

# Array to hold battery capacity in solutions
BatteryCapacityResults = zeros(2, nTurbineOpts)

# Loop through each turbine option
for turbine=1:nTurbineOpts

    # Read wind turbine power curve data. Data is also from NREL's SAM
    TurbineData = CSV.read(TurbineFiles[turbine])
    # Identify rated turbine power [kW]
    TurbinekW = maximum(TurbineData[2])

    # Wind cost data, based on scaling model from VanderMeer, Mueller-Stoffels, and Whitney 2017.
    # Overall numbers have been scaled to match Lazard
    # Upfront capital cost [$/kW]
    WindCapitalCost = exp(12.9 - 0.77 * log(TurbinekW) + 0.028 * log(TurbinekW)^2)/3.6
    # Annual cost [$/kw/year]
    WindAnnualCost = 35.00
    # Lifetime of wind project [years]
    WindLife = 20
    # Capital recovery factor
    CRFwind = LoanInterest * (1+LoanInterest)^WindLife / ((1+LoanInterest)^WindLife - 1)

    # Find generation of a single turbine at each time
    WindGeneration = zeros(nTimes)
    for t=1:nTimes
        if Speed[t] >= maximum(TurbineData[1])
            # If wind is too fast, the turbine cannot generate anything
            WindGeneration[t] = 0
        elseif Speed[t] <= minimum(TurbineData[1])
            # If wind is too slow, the turbine cannot generate anything
            WindGeneration[t] = 0
        else
            # Find generation at the current wind speed, interpolating between data points
            i = findlast(x -> x <= Speed[t], TurbineData[1])
            vi = TurbineData[1][i]
            vi1 = TurbineData[1][i+1]
            Pi = TurbineData[2][i]
            Pi1 = TurbineData[2][i+1]
            WindGeneration[t] = Pi + (Speed[t] - vi)*(Pi - Pi1)/(vi - vi1)
        end
    end

    for s=1:2

        # Try both fixed-tilt and single-axis tracking, and choose the appropriate
        # cost and generation parameters

        # Also print current parameters being used

        if s==1
            SolarGeneration = FixedGeneration
            SolarCapitalCost = FixedCapitalCost
            SolarAnnualCost = FixedAnnualCost
            println(TurbineOpts[turbine]," wind turbines, fixed-tilt solar")
        else
            SolarGeneration = TrackingGeneration
            SolarCapitalCost = FixedCapitalCost * TrackingCharge
            SolarAnnualCost = FixedAnnualCost * TrackingCharge
            println(TurbineOpts[turbine]," wind turbines, single-axis tracking solar")
        end

        #####################################
        ######### Define model ##############
        #####################################

        # Initialize solver
        m = Model(solver=CbcSolver())

        ###########################################
        ######### Decision Variables ##############
        ###########################################

        # Amount of solar capacity installed [kW]
        @variable(m, kWSolar >= 0)

        # Number of wind turbines installed (must be a whole number of turbines)
        @variable(m, nTurbines >= 0, Int)

        # Energy stored during each time step [kWh]
        @variable(m, StoreEnergyIn[1:nTimes] >= 0)

        # Energy removed from the battery during each time step [kWh]
        @variable(m, StoreEnergyOut[1:nTimes] >= 0)

        # The amount of energy stored in the battery at each time step [kWh]
        @variable(m, InStorage[1:nTimes+1] >= 0)

        # The storage of the battery [kWh]
        @variable(m, BatteryCapacity >= 0)

        ##############################################
        ######### Intermediate expressions ###########
        ##############################################

        # The total annualized cost of the chosen amount of solar [$/year]
        @expression(m, TotalSolarCost, kWSolar * (SolarCapitalCost * CRFsolar + SolarAnnualCost))

        # The total annualized cost of the chosen amount of wind [$/year]
        @expression(m, TotalWindCost, nTurbines * TurbinekW * (WindCapitalCost * CRFwind + WindAnnualCost))

        # The total battery cost. This includes the upfront cost of the battery,
        # as well as the annual cost of the battery degradation [$/yr]
        @expression(m, TotalBatteryCost, CostCapacity * BatteryCapacity * CRFbattery+sum(CostCapacity*StoreEnergyOut[t]/CycleLife for t = 1:nTimes))

        # The total power generated at a given time [kW]
        @expression(m, TotalGeneration[t=1:nTimes], WindGeneration[t] * nTurbines + SolarGeneration[t] * kWSolar)

        ###########################################
        ######### Objective function ##############
        ###########################################

        # Minimize the total cost of providing power for this grid [$/yr]
        @objective(m, Min, TotalSolarCost + TotalWindCost + TotalBatteryCost)

        ####################################
        ######### Constraints ##############
        ####################################

        # There must always be adequate power, and storage is limited to leftover power [kW]
        @constraint(m, [t=1:nTimes], TotalGeneration[t] + (1-RoundTripLosses) * StoreEnergyOut[t] - StoreEnergyIn[t]>= Load[t])

        # There must always be some spare capacity in case of sudden changes
        @constraint(m, [t=1:nTimes], TotalGeneration[t] + (1-RoundTripLosses) * InStorage[t] >= (1+PRM)*Load[t])

        # Storage conversion of energy constraint [kWh]
        @constraint(m, [t=1:nTimes], InStorage[t+1] == InStorage[t] + StoreEnergyIn[t] - StoreEnergyOut[t])

        # Maximum capacity constraint [kWh]
        @constraint(m, [t=1:nTimes], InStorage[t] <= BatteryCapacity)

        # Storage initialization constraint. While data are just for one year, it is
        # assumed that the cycle repeats afterwards, so the starting condition must
        # match the ending condition
        @constraint(m, InStorage[1] == InStorage[nTimes + 1])

        ########################################
        ######### Solve the model ##############
        ########################################

        solve(m)

        ###########################################
        ######### Display and save results ########
        ###########################################

        # Display values of key decision variables
        println("Number of wind turbines: ",getvalue(nTurbines))
        println("Solar kW: ", getvalue(kWSolar))
        println("Battery Capacity [kWh]: ", getvalue(BatteryCapacity))

        # Display cost, both as total annual cost and LCOE
        println("Total cost [\$/yr]: ", getobjectivevalue(m))
        println("LCOE [\$/kWh]: ", getobjectivevalue(m)/sum(Load))

        println("--------------------------")

        # Save results
        results = DataFrame(Month=SolarData[2],Day=Day,Time=SolarData[4],load=Load,Solar=(SolarGeneration*getvalue(kWSolar)),Wind=WindGeneration*getvalue(nTurbines),Stored=getvalue(InStorage)[1:nTimes],BattIn=getvalue(StoreEnergyIn),BattOut=getvalue(StoreEnergyOut))
        ResultsDF[s,turbine] = results

        CostResults[s, turbine] = getobjectivevalue(m)
        kWSolarResults[s, turbine] = getvalue(kWSolar)
        nTurbinesResults[s, turbine] = getvalue(nTurbines)
        BatteryCapacityResults[s, turbine] = getvalue(BatteryCapacity)

    end
end

println("--------------------------")
# Determine overall optimal results
OptimalPrice = minimum(CostResults)
OptimalCombination = findfirst(isequal(OptimalPrice), CostResults)
OptimalkWSolar = kWSolarResults[OptimalCombination]
OptimalnTurbines = nTurbinesResults[OptimalCombination]
OptimalBatteryCapacity = BatteryCapacityResults[OptimalCombination]

# Print results
println("Overall optimum:")
println("Turbine size: ",TurbineOpts[OptimalCombination[2]])
if OptimalCombination[1] == 1
    println("Fixed-tilt solar")
    println("PV Tilt Angle: ",rad2deg(CollectorTilt), "°")
    println("PV Azimuth Angle: ", rad2deg(CollectorAzimuth), "°")
else
    println("Single-axis tracking solar")
end
println("Number of wind turbines: ",OptimalnTurbines)
println("Solar kW: ", OptimalkWSolar)
println("Battery Capacity [kWh]: ", OptimalBatteryCapacity)

# Display cost, both as total annual cost and LCOE
println("Total cost [\$/yr]: ", OptimalPrice)
println("LCOE [\$/kWh]: ", OptimalPrice/sum(Load))

# Save hourly profile
CSV.write("results.csv",ResultsDF[OptimalCombination])
