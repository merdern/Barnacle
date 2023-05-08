
import pandas as pd
import numpy as np
import time
from scipy import signal
import streamlit as st
import plotly.graph_objs as go
import random


fs = 16

def DataProcessing(BarnData, ZerosData, low, high):

    
    
    #Define Calibration Coefficients, set for all Barnacles
    LDyn_0 = 1381
    Dynfit = [3.5734e-07, 5.1013e-06, 4.1913e-06, -0.0016, 0.9906]
    #LDyn = np.array([[-30,-20,-15,-10,-10,-5,0,0,0,0,5,10,10,15,20,30], [1.5648,1.2083,1.1237,1.0591,1.0178,0.9881,1,0.9721,1,1.0106,1.0018,1.0236,1.0359,1.0720,1.2350,1.7480]])
    Yawfit = [5.6971e-05, 2.2509e-04, 0.0531, -0.0661]
    #LYaw = np.array([[-30,-20,-15,-10,-10,-5,0,0,0,0,0,5,10,10,15,20,30],[-2.9570, -1.5197, -1.1453, -0.5479, -0.6622, -0.2944, -0.0473, -0.0340, -0.0402, -0.0449, -0.0666, 0.2183, 0.4996, 0.5528, 1.0175, 1.4180, 3.3176]])

    #Define other variables
    rho = 997
    vectstart = 795

    #Calibration curves
    yawcal1 = np.linspace(-45,45,91)
    yawcal2 = np.polyval(Yawfit, yawcal1)
    dyncal = np.polyval(Dynfit, yawcal1)*LDyn_0

    #Finding actual zeros data and removing from Barnacle Data
    zerosarray = np.array(ZerosData[low: high])
    zerosmean = np.mean(zerosarray, axis = 0)
    BarnData = BarnData - zerosmean[0:4]

    

    #Calibration interpolations
    denominator = np.mean(BarnData.iloc[:, 0:4], axis = 1)
    probeyaw = (BarnData.iloc[:,1]-BarnData.iloc[:,3])/denominator
    probepitch = (BarnData.iloc[:,0]-BarnData.iloc[:,2])/denominator
    yaw = np.interp(probeyaw, yawcal2, yawcal1, period = 360)
    pitch = np.interp(probepitch, yawcal2, yawcal1, period = 360)

    #Dynamic calculations
    pitchbigger = (abs(probepitch)>abs(probeyaw))*1
    probemax = pitchbigger*probepitch + (pitchbigger -1)*probeyaw
    dyn = np.interp(probemax, yawcal1, dyncal)

    #Find and resolve velocities
    U = np.sqrt(2 *-dyn * np.mean(BarnData.iloc[:, 0:4], axis = 1)/rho)
    U = np.real(U)

    #Removing noise from data, not sure if this works with live data
    #sos = signal.butter(8,0.5,output ='sos')
    #U = signal.sosfiltfilt(sos, U)
    #plt.plot(U)
    
    #Resolve velocities in x, y and z directions
    Ux = U*np.cos(pitch*np.pi/180)*np.cos(yaw*np.pi/180)
    Uy = U*np.cos(pitch*np.pi/180)*np.sin(yaw*np.pi/180)
    Uz = U*np.sin(pitch*np.pi/180)
    
    

    return(yaw,pitch,U,Ux,Uy,Uz)


def findintensity(U, Ux, Uy, Uz):

         
    #Find Turbulence Intensity in all dimensions
    intensity = round((np.std(U))/(np.mean(U))*100,2)
    intensityx = round((np.std(Ux))/(np.mean(U))*100,2)
    intensityy = round((np.std(Uy))/(np.mean(U))*100,2)
    intensityz = round((np.std(Uz))/(np.mean(U))*100,2)
        
    return(intensity, intensityx, intensityy, intensityz)

def findlengthscale(U):

    #Use an Autocorrelation to estimate a length scale for the data
    Up = U-np.mean(U)
    please = np.correlate(Up,Up, mode='full')
    please = please[round(len(please)/2):]
    asign = np.sign(please)
    signchange = ((np.roll(asign, 1) - asign) !=0). astype(int)
    signchange[0] = 0
    id = np.nonzero(signchange)  
    id = int(id[0][0])
    please = please[0:id]
    integral = np.trapz(please)/fs/len(please) * np.mean(U)

    return(integral)


def godata():
    #Data! was clicked
    st.session_state.godata = not st.session_state.godata
    st.session_state.data = Barndata

def stopdata():
     st.session_state.godata = False
     st.session_state.data = Barndata

def callback():
    #Button was clicked
    st.session_state.button_clicked = not st.session_state.button_clicked
    st.session_state.low = low
    st.session_state.high = high
    st.session_state.zerosdata = ZerosData

def refreshdata():
    st.session_state.data = pd.DataFrame(columns = ['col1','col2','col3','col4',])
    st.session_state.godata = False
    

#Define page layout and title
st.set_page_config(page_title='Barnacle Dashboard', page_icon='🌊', layout='wide', initial_sidebar_state= 'expanded')
st.title('Barnacle Dashboard 📈')

#Upload and read a Zero flow file
ZerosFile = st.sidebar.file_uploader('Upload a file with zero flow data')

placeholderbutton = st.sidebar.empty()
placeholderdata = st.sidebar.empty()
placeholder = st.empty()

chart_layout = go.Layout(
                    title='Velocity over Time',
                    xaxis=dict(title='Time'),
                    yaxis=dict(title='Velocity'),
                    hovermode='closest',
                    legend = dict(x = 0, y = 1, font=dict(size=16)),
                    height = 700
                )

chart_layoutpitch = go.Layout(
                    title='Pitch over Time',
                    xaxis=dict(title='Time'),
                    yaxis=dict(title='Pitch Angle'),
                    hovermode='closest',
                    height = 400
                )

chart_layoutyaw = go.Layout(
                    title='Yaw over Time',
                    xaxis=dict(title='Time'),
                    yaxis=dict(title='Yaw Angle'),
                    hovermode='closest',
                    height = 400
                )
toplot = st.sidebar.multiselect('Select Velocity Component(s) to plot on axis', options = ['Velocity Magnitude', 'Ux', 'Uy', 'Uz'])

if "godata" not in st.session_state:
    st.session_state.godata = False

if "traceU" not in st.session_state:
    st.session_state.traceU = []
if "traceux" not in st.session_state:
    st.session_state.traceux = []
if "traceuy" not in st.session_state:
    st.session_state.traceuy = []
if "traceuz" not in st.session_state:
    st.session_state.traceuz = []
if "figpitch" not in st.session_state:
    st.session_state.figpitch = []
if "figyaw" not in st.session_state:
    st.session_state.figyaw = []
if "length_scale" not in st.session_state:
    st.session_state.length_scale = []

low, high = [], []
if ZerosFile is not None:

    ZerosData = pd.read_table(ZerosFile, sep = ',')
    if "button_clicked" not in st.session_state:
            st.session_state.button_clicked = False
            st.session_state.zerosdata = []
    
    with placeholderbutton.container():
        if st.session_state.button_clicked == False:
                low, high = st.slider('Select range of Slack Water:', 0, len(ZerosData), (1300, len(ZerosData)))
                st.button('Confirm the Slack Water Time are correct', on_click=callback, use_container_width= True)
                st.write('Upload Zero Flow Data and initialise data stream using the sidebar')            

        else:
                st.button('Set Zero Flow region', on_click= callback,  use_container_width= True)

    if "data" not in st.session_state:
        Barndata = pd.DataFrame(columns = ['col1','col2','col3','col4',])
    else: Barndata = st.session_state.data

    
    if "godata" not in st.session_state:
         st.session_state.godata = False            

    with placeholderdata.container():
        col1, col2 = st.columns(2)
        with col1:
            st.button('Data!', on_click = godata, use_container_width= True)
        with col2:
             st.button('Stop Data', on_click = stopdata, use_container_width=True)
             
        st.button('Refresh Data', on_click = refreshdata, use_container_width= True)

        
    

    Barndata = [random.randint(350, 400) for _ in range(4)]
    Barndata = [Barndata/100 for Barndata in Barndata]
    Barndata = pd.DataFrame([Barndata])


    while st.session_state.godata:
                # Simulate data stream
                newdata = np.random.rand(1, 4) - 0.5
                newdata = pd.DataFrame(newdata)
                time.sleep(1/16)
                if 4.5 > (Barndata.iloc[-1] + newdata.iloc[0]).any() > 3.6:
                    Barndata.loc[len(Barndata)] = Barndata.iloc[-1] + newdata.iloc[0]
                else:
                     Barndata.loc[len(Barndata)] = Barndata.iloc[-1] - newdata.iloc[0]           
                
                
               

                if len(Barndata) != 0 and len(Barndata) % 16 == 0:

                    yaw,pitch,U,Ux,Uy,Uz = DataProcessing(Barndata, st.session_state.zerosdata, st.session_state.low, st.session_state.high)
                    times = np.arange(len(Ux)) / 16
                    intensity, intensityx, intensityy, intensityz = findintensity(U, Ux, Uy, Uz)
                    length_scale = findlengthscale(U)
                    length_scale = np.abs(round(length_scale, 2))
                    st.session_state.length_scale = length_scale

                    trace_U = go.Scatter(x = times, y = U, name = f"U - Intensity = {intensity}%")
                    trace_ux = go.Scatter(x = times, y = Ux, name = f"Ux - Intensity = {intensityx}%")
                    trace_uy = go.Scatter(x = times, y = Uy, name = f"Uy - Intensity = {intensityy}%")
                    trace_uz = go.Scatter(x = times, y = Uz, name = f"Uz - Intensity = {intensityz}%")
                    st.session_state.traceU = trace_U
                    st.session_state.traceux = trace_ux
                    st.session_state.traceuy = trace_uy
                    st.session_state.traceuz = trace_uz
                    
                    trace_pitch = go.Scatter(x = times, y=pitch, name='Pitch', mode='lines')
                    fig_pitch = go.Figure(data = trace_pitch, layout=chart_layoutpitch)
                    st.session_state.figpitch = fig_pitch
                                                
                    trace_yaw = go.Scatter(x = times, y=yaw, name='Yaw', mode='lines')
                    fig_yaw = go.Figure(data = trace_yaw, layout=chart_layoutyaw)
                    st.session_state.figyaw = fig_yaw

                    with placeholder.container():

                        col1, col2 = st.columns(2, gap = 'medium')

                        with col1:

                                st.subheader(f"Estimated length scale for this data is {length_scale}m")

                                                               
                                traces_vel = []
                                for option in toplot:
                                    if option == 'Velocity Magnitude' :
                                        traces_vel.append(st.session_state.traceU)
                                    elif option == 'Ux':
                                        traces_vel.append(st.session_state.traceux)
                                    elif option == 'Uy':
                                        traces_vel.append(st.session_state.traceuy)
                                    elif option == 'Uz':
                                        traces_vel.append(st.session_state.traceuz)
                                if len(traces_vel) != 0:
                                    fig_vel = go.Figure(data=traces_vel, layout=chart_layout)
                                    fig_vel.update_layout(showlegend = True)
                                    st.plotly_chart(fig_vel, use_container_width=True)
                        
                        with col2:
                            st.plotly_chart(fig_pitch, use_container_width= True)
                            st.plotly_chart(fig_yaw,  use_container_width= True)

    with placeholder.container():

                        col1, col2 = st.columns(2, gap = 'medium')

                        with col1:

                                st.subheader(f"Estimated length scale for this data is {st.session_state.length_scale}m")

                                                 
                                traces_vel = []
                                for option in toplot:
                                    if option == 'Velocity Magnitude' :
                                        traces_vel.append(st.session_state.traceU)
                                    elif option == 'Ux':
                                        traces_vel.append(st.session_state.traceux)
                                    elif option == 'Uy':
                                        traces_vel.append(st.session_state.traceuy)
                                    elif option == 'Uz':
                                        traces_vel.append(st.session_state.traceuz)
                                
                                fig_vel = go.Figure(data=traces_vel, layout=chart_layout)
                                fig_vel.update_layout(showlegend = True)
                                st.plotly_chart(fig_vel, use_container_width=True)
                        
                        with col2:
                            st.plotly_chart(st.session_state.figpitch, use_container_width= True)
                            st.plotly_chart(st.session_state.figyaw,  use_container_width= True)
                        
