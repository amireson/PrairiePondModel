# Pond water-isotope mass balance model: Plotting script
# 
# Edward Bam(1,2)<br>
# Andrew Ireson(1,3)
# 
# (1) Global Institute for Water Security, University of Saskatchewan www.usask.ca/GIWS <br>
# (2) e.bam@usask.ca <br>
# (3) andrew.ireson@usask.ca   http://usask.ca/~andrew.ireson
# 

from matplotlib import pyplot as plt
import numpy as np
import pandas
import matplotlib.dates as mdates

# Load observations (heavy)
print "Loading data\n"
pond=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Isotopes')
depth=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Depth')
pardata = pandas.read_excel('Input.xlsx',sheet_name='Pars',index_col='Parameter')
df=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Input')
output=pandas.read_excel('Output.xlsx',index_col=0,sheet_name='Timeseries')

print "Processing plot\n"

StartDate=pardata['Value'].loc['StartDate']
Year=StartDate.year

# Get Plot parameters
PlotStartDate=pardata['Value'].loc['PlotStart']
PlotEndDate=pardata['Value'].loc['PlotEnd']
TimeStepPlot=pardata['Value'].loc['TimeStep']
sp1_bot=pardata['Value'].loc['sp1_bot']
sp1_top=pardata['Value'].loc['sp1_top']
sp2_bot=pardata['Value'].loc['sp2_bot']
sp2_top=pardata['Value'].loc['sp2_top']
sp3_bot=pardata['Value'].loc['sp3_bot']
sp3_top=pardata['Value'].loc['sp3_top']
sp4_bot=pardata['Value'].loc['sp4_bot']
sp4_top=pardata['Value'].loc['sp4_top']
sp1_del=pardata['Value'].loc['sp1_del']
sp2_del=pardata['Value'].loc['sp2_del']
sp3_del=pardata['Value'].loc['sp3_del']
sp4_del=pardata['Value'].loc['sp4_del']

# Get time step in days
df['dt']=0.
for i in range(1,len(df.index)):
    df.ix[i,'dt']=(df.index[i]-df.index[i-1]).days

# Convert precipitation into mm/d
df['P'] = df['P']/df['dt']

dfPrecip=df.copy()

S=PlotStartDate
E=PlotEndDate

tS=output.index[0]
tE=output.index[-1]
depthP=depth[tS:tE]

xt=pandas.date_range(S,E,freq=TimeStepPlot)
xtl=xt.strftime('%d-%b')

f, (ax0,ax1,ax2,ax3)=plt.subplots(4, sharex=True, figsize=(7,6))

# Precipitation plot:
t=dfPrecip.index
P=(dfPrecip['P']).values
ax0.step(t,P,color='k',label='Precipitation intensity')
ax0.set_ylabel('mm/d',fontsize=12)
ax0.legend(loc=1,fontsize=12)
ax0.set_ylim((sp1_bot,sp1_top))
ax0.set_xlim(S,E)
yt=np.arange(sp1_bot,sp1_top+sp1_del,sp1_del)
ax0.set_yticks(yt)
ax0.grid()
ax0.get_xaxis().set_tick_params(direction='in')
  
# Pond water level
t1=depth.index
d1=depth['Depth'].values
t2=depthP.index
d2=depthP['Depth'].values
t3=output.index
d3=output['k'].values

if len(depthP['Depth'])>30:
    ms='.'
else:
    ms='o'

ax1.plot(t1,d1,ms,color='b',label='Observations')
ax1.plot(t3,d3,'-r',label='Simulation',linewidth=2.5)

ax1.get_xaxis().set_tick_params(direction='in')

ax1.set_ylim((sp2_bot,sp2_top))
yt=np.arange(sp2_bot,sp2_top,sp2_del)
ax1.set_yticks(yt)
ax1.grid()
ax1.legend(loc=3,fontsize=12,bbox_to_anchor=(0, -0.05))

ax1.set_ylabel('Water level (m)',fontsize=12)

# O18 plot
tO1=pond['O18'].index
O1=pond['O18'].values
tO2=output['O18'].index
O2=output['O18'].values
ax2.plot(tO1,O1,'o',color='b')
ax2.plot(tO2,O2,'-r',linewidth=2.5)
ax2.set_ylim((sp3_bot,sp3_top))
yt=np.arange(sp3_bot,sp3_top,sp3_del)
ax2.set_yticks(yt)
ax2.grid()
ax2.get_xaxis().set_tick_params(direction='in')
ax2.set_ylabel(u'$\delta ^{18} O  (\u2030) $',fontsize=12)

# H2 plot
tH1=pond['H2'].index
H1=pond['H2'].values
tH2=output['H2'].index
H2=output['H2'].values
ax3.plot(tH1,H1,'o',color='b')
ax3.plot(tH2,H2,'-r',linewidth=2.5)
ax3.set_ylim((sp4_bot,sp4_top))
ax3.grid()
ax3.set_ylabel(u'$\delta ^{2} H  (\u2030) $',fontsize=12)
yt=np.arange(sp4_bot,sp4_top,sp4_del)
ax3.set_yticks(yt)

ax3.set_xticks(xt)  
ax3.set_xticklabels(xtl)
ax3.set_xlabel('Date (%d)' % Year,fontsize=12)


ax0.set_position([0.12, 0.85, 0.85, 0.1])
ax1.set_position([0.12, 0.6, 0.85, 0.24])
ax2.set_position([0.12, 0.35, 0.85, 0.24])
ax3.set_position([0.12, 0.1, 0.85, 0.24])

# Save figure file
plt.savefig('Figure.png',dpi=300)

