# import libraries 
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.dates as mdates 

#Monthly Data
case1m = pd.read_csv('pvwatts_monthlyAug5.csv')
case2m = pd.read_csv('pvwatts_monthlyFeb9.csv')
case3m = pd.read_csv('pvwatts_monthlyWest.csv')
#Hourly Data
case1h = pd.read_csv('pvwatts_hourlyAug5.csv') #Aug5 angle
case2h = pd.read_csv('pvwatts_hourlyFeb9.csv') #Feb9 angle
case3h = pd.read_csv('pvwatts_hourlyWest.csv') #West facing Aug5 angle

#Hourly AC System Output for Aug 5-- need month=8 day=5 hour=0:23
plt.figure(figsize=(9, 5))

# Filter data for month=8, day=5, and hour=0:23
#df1 is for Case 1 specifically; df2 for Case 2, etc!
df1 = case1h[(case1h['Month'] == 2) &  
                 (case1h['Day'] == 9) &    
                 (case1h['Hour'] >= 0) & (case1h['Hour'] <= 23)]   # Hour range 0-23

df2 = case2h[(case2h['Month'] == 2) &  
                 (case2h['Day'] == 9) &    
                 (case2h['Hour'] >= 0) & (case2h['Hour'] <= 23)]   # Hour range 0-23

df3 = case3h[(case3h['Month'] == 2) &  
                 (case3h['Day'] == 9) &    
                 (case3h['Hour'] >= 0) & (case3h['Hour'] <= 23)]   # Hour range 0-23

# Plot the filtered data
plt.plot(df1['Hour'], df1['AC System Output (W)'], label='Case 1')
plt.plot(df2['Hour'], df2['AC System Output (W)'], label='Case 2')
plt.plot(df3['Hour'], df3['AC System Output (W)'], label='Case 3')

#plt.plot(case1h['Hour'], case1h['AC System Output (W)'], label='Case 1')
#plt.plot(case2m['Month'], case2m['Daily Average POA Irradiance (kWh/m2/day)'], label='Case 2')
#plt.plot(case3m['Month'], case3m['Daily Average POA Irradiance (kWh/m2/day)'], label='Case 3')

plt.title('Hourly AC System Output on February 9')
plt.xlabel('Hour')
plt.ylabel('AC System Output (W)')
plt.xticks(case1h['Hour'])

# Show the plot
plt.grid(True)
plt.legend()
#plt.tight_layout()
plt.show()