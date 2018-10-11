# # Pond water-isotope mass balance model
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
from scipy import optimize

# Load observations (heavy)
print "Loading observed data\n"
pond=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Isotopes')
depth=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Depth')

# Function to read input file:
def ReadInput(E,O,H,T,A_H2,A_O18,StartDate):
    # Use a library/dataframe for time series data:
    # Read input data from Template
    df=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Input')

    # Get time step in days
    df['dt']=0.
    for i in range(1,len(df.index)):
        df['dt'][i]=(df.index[i]-df.index[i-1]).days
        
    # Convert precipitation into mm/d
    df['P'] = df['P']/df['dt']

    # Convert relative humidity into a fraction
    df['H']=H/100.
    df['T']=T

    # Initialize columns to be populated:
    df['k']=0.

    # Create dataframes for the isotopes
    df_O18=pandas.DataFrame()
    df_O18['Prec']=df['P_O18'].copy()
    df_O18['Atm']=A_O18
    df_O18['Pond']=0.

    df_H2=pandas.DataFrame()
    df_H2['Prec']=df['P_H2'].copy()
    df_H2['Atm']=A_H2
    df_H2['Pond']=0.
    
    return df.loc[df.index>StartDate], df_O18.loc[df.index>StartDate], df_H2.loc[df.index>StartDate]

# Core model functions
def PondDAV(pars,k=0.,V=0.,A=0.):
    beta=pars['beta']
    s=pars['s']
    ki=pars['ki']
    if k!=0:
        A=s*(k/ki)**(2/beta)
        V=A*k/(1+2/beta)
    elif A!=0:
        k=ki*(A/s)**(beta/2)
        V=A*k/(1+2/beta)
    elif V!=0:
        k=(V*(1+2/beta)/s*ki**(-beta/2))**(1/(1+2/beta))
        A=s*(k/ki)**(2/beta)
    return k,A,V
PondDAV=np.vectorize(PondDAV)

# Fractionation factor functions
def Alpha_O18(T):
    T=T+273
    Al = -7.685+6.7123e3/T-1.6664e6/T**2+0.35041e9/T**3
    Al = np.exp(Al/1000)
    return Al

def Alpha_H2(T):
    T=T+273
    Al = 1158.8e-9*T**3 - 1620.1e-6*T**2 + 794.84e-3*T - 161.04 + 2.9992e9/T**3
    Al = np.exp(Al/1000)
    return Al

# Isotope balance model:
# t = time (d)
# dt = current time step length (d)
# df['P'] = Precip during timestep (mm/d)
# df['E'] = Evap during timestep(mm/d)
# df['O'] = Pond infiltration during timestep (mm/d)
# df['T'] = Temperature during time step (deg C)
# df['del_i'] = Initial isotope conc. of pond
# df['del_a'] = atmospheric isotopic conc. during timestep
# df['H'] = Relative humidity during timestep (-)
# df['k'] = Water level.
# pars = library of parameters, including:
    # s = Pond geom. par
    # beta = Pond geom. par
    # ki = Pond geom. par
    # thetanC = coefficient for kinetic isotopic separation

def ModelFun(E,O,k_i,del_i,dt,df,df_iso,pars,p_iso,alphafun):
    
    # Water balance
    k_t=k_i+(df['P']-E-O)/1000.*dt
    dum,A,V=PondDAV(pars,k_i)
    
    # Isotope balance
    Al=alphafun(df['T'])
    Ep=1000*(Al-1)
    DelEp=p_iso['thetanC']*(1-df['H'])
    
    B=(df['H']*df_iso['Atm']+DelEp/1.+Ep/Al)/(1-df['H']+DelEp/1000.)
    #B=(df['H']*df_iso['Atm']+DelEp/1000.+Ep/Al)/(1-df['H']+DelEp/1000.)
    C=(df['H']-DelEp/1000.-Ep/1000./Al)/(1-df['H']+DelEp/1000.)
    
    X1=np.exp(-(A*df['P']/1000.+A*E/1000.*C)*dt/V)
    X2=(E/1000.*B+df['P']/1000.*df_iso['Prec'])/(df['P']/1000.+E/1000.*C)
    X3=(1-np.exp(-(A*df['P']/1000.+A*E/1000.*C)*dt/V))
    del_t=del_i*X1+X2*X3
    
    if k_t<0:
        k_t=0.
        del_t=0.

    return k_t,del_t

def RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2):

    # Solve model
    df['k'][0]=k_ini
    df_O18['Pond'][0]=O18_ini
    df_H2['Pond'][0]=H2_ini
    for i in range(1,len(df.index)):
        dt=(df.index[i]-df.index[i-1]).days
        k_t,del_O18_t = ModelFun(E,O,df['k'][i-1],df_O18['Pond'][i-1],df['dt'][i],df.iloc[i],df_O18.iloc[i],pars,p_O18,Alpha_O18)
        df['k'][i]=k_t
        df_O18['Pond'][i]=del_O18_t
        k_t,del_H2_t = ModelFun(E,O,df['k'][i-1],df_H2['Pond'][i-1],df['dt'][i],df.iloc[i],df_H2.iloc[i],pars,p_H2,Alpha_H2)
        df_H2['Pond'][i]=del_H2_t

    dum,df['A'],df['V']=PondDAV(pars,k=df['k'])
    
    output=df['k'].to_frame()
    output['O18']=df_O18['Pond'].copy()
    output['H2']=df_H2['Pond'].copy()
    
    return output, df, df_O18, df_H2


# Functions for stage one - optimize the total water loss to the water levels
def OptFun_Level(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2):
    output, df, df_O18, df_H2 = RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2)
    RMSE=ObFun_Level(depth,output)
    #print Popt, RMSE
    return RMSE

def ObFun_Level(depth,output):
    t_o=depth.index
    k_o=depth.values.squeeze()
    t_s=output.index
    k_s=output['k'].values.squeeze()

    if len(t_o) > len(t_s):
        t=t_s.to_julian_date()
        tx=t_o.to_julian_date()
        tx=tx-t[0]
        t=t-t[0]
        k_o=np.interp(t,tx,k_o)
    else:
        t=t_o.to_julian_date()
        tx=t_s.to_julian_date()
        tx=tx-t[0]
        t=t-t[0]
        k_s=np.interp(t,tx,k_s)

    RMSE = np.sqrt(np.mean((k_o-k_s)**2.))
    return RMSE

# Functions for stage two - optimize the evaporation to isotope values
def OptFun_Iso(E,EI,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2):
    O = EI - E
    
    output, df, df_O18, df_H2 = RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2)
    RMSE_O18=ObFun_Iso(pond,output,'O18')/np.abs(pond['O18'].mean())
    RMSE_H2=ObFun_Iso(pond,output,'H2')/np.abs(pond['H2'].mean())
    RMSE=RMSE_O18+RMSE_H2
    return RMSE

# RMSE for the pond depth function:
def ObFun_Iso(pond,output,isostr):
    t_o=pond.index
    iso_o=pond[isostr].values.squeeze()
    t_s=output.index
    iso_s=output[isostr].values.squeeze()
    
    if len(t_o) > len(t_s):
        t=t_s.to_julian_date()
        tx=t_o.to_julian_date()
        tx=tx-t[0]
        t=t-t[0]
        iso_o=np.interp(t,tx,iso_o)
    else:
        t=t_o.to_julian_date()
        tx=t_s.to_julian_date()
        tx=tx-t[0]
        t=t-t[0]
        iso_s=np.interp(t,tx,iso_s)


    RMSE = np.sqrt(np.mean((iso_o-iso_s)**2.))
    return RMSE

###############################################################################################################
# Function to plot parameter sensitivity
def SensiContourPlot(fn,levels,labels,mytitle):
    f=open(fn,'r')
    P=[]
    k=np.zeros((9,3))
    Z=np.zeros((9,3))
    for i in range(3):
        P.append(float(f.readline().strip()))
        for j in range(9):
            myline=f.readline().strip().split()
            k[j,i]=float(myline[1])
            Z[j,i]=float(myline[2])
    EI=np.array([0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9])

    cmap=plt.cm.get_cmap('hot')
    xt=[0.6, 1.0, 1.4]
    fs=16

    fig=plt.figure(figsize=(5,6))
    ax = fig.add_axes([0.15, 0.23, 0.84, 0.7])
    x=np.array([0.4, 0.8])
    z=np.array([Z[:,0],Z[:,0]]).transpose()
    plt.contourf(x,EI,z,levels,cmap=cmap)
    x=np.array([0.8, 1.2])
    z=np.array([Z[:,1],Z[:,1]]).transpose()
    plt.contourf(x,EI,z,levels,cmap=cmap)
    x=np.array([1.2, 1.6])
    z=np.array([Z[:,2],Z[:,2]]).transpose()
    plt.contourf(x,EI,z,levels,cmap=cmap)

    ax.set_xticks(xt)
    ax.set_xticklabels(labels,fontsize=fs)
    ax.set_yticklabels(EI,fontsize=fs)
    plt.plot([0.8,0.8],[0.1,0.9],'-k')
    plt.plot([1.2,1.2],[0.1,0.9],'-k')
    plt.title(mytitle,fontsize=fs)
    plt.ylabel('E/(E+O)',fontsize=fs)

    cbaxes = fig.add_axes([0.1, 0.1, 0.8, 0.04]) 
    # cb = plt.colorbar(ax1, cax = cbaxes)  
    cbar=plt.colorbar(cax=cbaxes,orientation='horizontal')
    # cbar.ax.tick_params(labelsize=fs)
    cbar.set_label(label='del O18',fontsize=fs)
    
    plt.savefig(fn.replace('txt','png'),dpi=300)


# Read parameters
pardata = pandas.read_excel('Input.xlsx',sheet_name='Pars',index_col='Parameter')

StartDate=pardata['Value'].loc['StartDate']
E=pardata['Value'].loc['E']
O=pardata['Value'].loc['O']

# Pond geometry parameters:
pars={}
pars['s']=pardata['Value'].loc['s']
pars['beta']=pardata['Value'].loc['beta']
pars['ki']=pardata['Value'].loc['ki']

ShallIOptimize=pardata['Value'].loc['Optimize']

# Isotopic parameters:
p_O18={}
p_H2={}
p_O18['thetanC']=pardata['Value'].loc['O18_thetanC']
p_H2['thetanC']=pardata['Value'].loc['H2_thetanC']

# Initial conditions:
k_ini=np.interp(StartDate.to_julian_date(),depth.index.to_julian_date(),depth['Depth'].values)
O18_ini=np.interp(StartDate.to_julian_date(),pond.index.to_julian_date(),pond['O18'].values)
H2_ini=np.interp(StartDate.to_julian_date(),pond.index.to_julian_date(),pond['H2'].values)
    
# Read temperature and humidity:
H=pardata['Value'].loc['H']
T=pardata['Value'].loc['T']
    
# Read atmospheric isotopic composition:
A_H2=pardata['Value'].loc['A_H2']
A_O18=pardata['Value'].loc['A_O18']

# Read timeseries input data
df,df_O18,df_H2=ReadInput(E,O,H,T,A_H2,A_O18,StartDate)

# Code block to run the model in one of three possible configurations:
if ShallIOptimize==0:
    # Run model with given parameters, print output to screen
    output, df, df_O18, df_H2 = RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2)
    print "    DRIVING DATA AND WATER BALANCE:"
    print df
    print ""
    print "    SIMULATED O18:"
    print df_O18
    print ""
    print "    SIMULATED H2:"
    print df_H2

elif ShallIOptimize==1:
    # Optimize model to pond level and pond isotopic concentrations
    print "Stage one optimization: get total water loss from the pond"
    EI = E + O
    O = 0.
    EI=optimize.fmin(func=OptFun_Level, x0=E,args=(O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2,))

    output, df, df_O18, df_H2 = RunModel(EI,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2)
    RMSE1=ObFun_Level(depth,output)

    print "Stage two optimization: get evaporative loss from the pond"
    E=optimize.fmin(func=OptFun_Iso, x0=E,args=(EI,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2,))
    #E = np.min((EI,np.max((0,E))))
    O = EI - E
    output, df, df_O18, df_H2 = RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,df,df_O18,df_H2,pars,p_O18,p_H2)
    RMSE_O18=ObFun_Iso(pond,output,'O18')/np.abs(pond['O18'].mean())
    RMSE_H2=ObFun_Iso(pond,output,'H2')/np.abs(pond['H2'].mean())
    RMSE2=RMSE_O18+RMSE_H2

    print "RMSE (level) = %0.2f" % RMSE1
    print "RMSE (%s) = %0.2f" % (ShallIOptimize,RMSE2)
    print "E = %0.2f" % E
    print "O = %0.2f" % O

    numday=df.index[-1].to_julian_date()-df.index[0].to_julian_date()
    k_fin=df['k'][-1]
    outEI={}
    outEI['RMSE_level']=RMSE1
    outEI['RMSE_iso']=RMSE2
    outEI['E (mm/d)']=E
    outEI['O (mm/d)']=O
    outEI['E (mm)']=E*numday
    outEI['O (mm)']=O*numday
    outEI['Number of days']=numday
    P=(df['P'][1:]*df['dt'][1:]).sum()/numday
    outEI['P (mm/d)']=P
    outEI['P (mm)']=P*numday
    outEI['Sim_deltaStorage (mm)']=(E+O-P)*numday
    outEI['Obs_deltaStorage (mm)']=1000*(k_ini-k_fin)
    outEI['k ini (mm)']=k_ini*1000
    outEI['k fin (mm)']=k_fin*1000
    outEI['Pond area ini (m2)']=PondDAV(pars,k=k_ini)[1]
    outEI['Pond area fin (m2)']=PondDAV(pars,k=k_fin)[1]
    outEI=pandas.DataFrame(data=outEI,index=np.array(['Simulation:']))


    # Save result
    writer=pandas.ExcelWriter('Output.xlsx')
    outEI.transpose().to_excel(writer,'Outputs')
    output.to_excel(writer,'Timeseries')
    writer.save()
   
elif ShallIOptimize==2:
    # Perform a simple sensitivity analysis
    
    sensi=pandas.read_excel('Input.xlsx',index_col=0,sheet_name='Sensi',convert_float=False)
    n=len(sensi)
    
    # Loop through each parameter to be perturbed:
    for i in range(n):
        # Set up a default parameter matrix with all baseline parameter values.
        middle=np.tile(sensi['Middle'].values,(3,1))
        ParMat=pandas.DataFrame(middle,index=['Lower','Middle','Upper'],columns=sensi['Parameter'].values)
        # Modify values of selected parameter
        ParName=sensi.iloc[i].Parameter
        ParMat[ParName]['Lower']=sensi['Lower'].iloc[i]
        ParMat[ParName]['Middle']=sensi['Middle'].iloc[i]
        ParMat[ParName]['Upper']=sensi['Upper'].iloc[i]
        
        fn=[str(sensi.index[i])+'Sensitivity.txt'][0]
        print fn
        print(ParMat)

        f=open(fn,'w')

        Label1=sensi.iloc[i].Label1.strip('"')
        Label2=sensi.iloc[i].Label2.strip('"')
        Label3=sensi.iloc[i].Label3.strip('"')
        labels=[Label1, Label2, Label3]
        mytitle=sensi.iloc[i].Title.strip('"')
        # Perform sensitivity analysis
        EI=E+O
        par_values=sensi[['Lower','Middle','Upper']].iloc[i].values
        for j in range(3):
            f.write('%.3f\n'%float(par_values[j]))
            for EvapFraction in np.linspace(0.1,0.9,9):
                # Assign parameters:
                Hfact=ParMat.iloc[j]['Hfact']
                precipfact=ParMat.iloc[j]['precipfact']
                Tfact=ParMat.iloc[j]['Tfact']
                dmfact=ParMat.iloc[j]['dmfact']
                betafact=ParMat.iloc[j]['betafact']
                sfact=ParMat.iloc[j]['sfact']
                
                # Run model
                E=EI*EvapFraction
                O=EI*(1.0-EvapFraction)
                
                dfM=df.copy()
                dfM_O18=df_O18.copy()
                dfM_H2=df_H2.copy()
                dfM['P']=dfM['P']*precipfact
                dfM['T']=dfM['T']*Tfact
                dfM['H']=dfM['H']*Hfact
                dfM_O18['Atm']=dfM_O18['Atm']*dmfact
                dfM_H2['Atm']=dfM_H2['Atm']*dmfact
                
                parsM=pars.copy()
                parsM['beta']=pars['beta']*betafact
                parsM['s']=pars['s']*sfact

                output, dfM, dfM_O18, dfM_H2 = RunModel(E,O,k_ini,H2_ini,O18_ini,StartDate,dfM,dfM_O18,dfM_H2,parsM,p_O18,p_H2)
                
                f.write('%.2f %.4f %.4f\n'%(EvapFraction, dfM['k'].values[-1], dfM_O18['Pond'].values[-1]))
        levels=np.linspace(-12,-4,9)
        f.close()
        SensiContourPlot(fn,levels,labels,mytitle)


