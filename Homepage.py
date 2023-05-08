

#TO CONVERT TO EXE 
#https://youtu.be/3wZ7GRbr91g


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import time
from scipy import signal
import streamlit as st
import plotly.graph_objs as go
import random
import io
from itertools import groupby
import datetime


fs = 16
#BarnFile = "C:\\Users\\twcy2P\\OneDrive\\Documents\\Work\\Uni\\Third Year\\GBPD\MATLAB\\Mon1501"
#ZerosFile = "C:\\Users\\twcy2\\OneDrive\\Documents\\Work\\Uni\\Third Year\\GBPD\\MATLAB\\Mon1527"
#BarnData = pd.read_table(BarnFile, sep = ',')
#ZerosData = pd.read_table(ZerosFile, sep = ',')


@st.cache_data
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

    sos = signal.butter(8, 0.5, output ='sos')
    BarnDataFilt = np.zeros_like(BarnData)
    BarnDataFilt = pd.DataFrame(index=BarnData.index, columns=BarnData.columns)
    for column in BarnData.columns:
        BarnDataFilt[column] = signal.sosfiltfilt(sos, BarnData[column])

    #Finding actual zeros data and removing from Barnacle Data
    zerosarray = np.array(ZerosData[low: high])
    zerosmean = np.mean(zerosarray, axis = 0)
    BarnData = BarnDataFilt - zerosmean

    #Calibration interpolations
    denominator = np.mean(BarnData.iloc[:, 0:4], axis = 1)
    probeyaw = (BarnData.iloc[:,1]-BarnData.iloc[:,3])/denominator
    probepitch = (BarnData.iloc[:,0]-BarnData.iloc[:,2])/denominator
    #yaw = np.interp(probeyaw, yawcal2, yawcal1, period = 360)
    #pitch = np.interp(probepitch, yawcal2, yawcal1, period = 360)

    # find the rows where both probeyaw and probepitch are within the calibration range
    in_range = ((probeyaw >= yawcal2[25]) & (probeyaw <= yawcal2[65]) &
    (probepitch >= yawcal2[25]) & (probepitch <= yawcal2[65]))

    dyn = np.full_like(probeyaw, np.nan)
    yaw = np.full_like(probeyaw, np.nan)
    pitch = np.full_like(probepitch, np.nan)
# apply the interpolation only for the in-range rows
    yaw[in_range] = np.interp(probeyaw[in_range], yawcal2, yawcal1, period=360)
    pitch[in_range] = np.interp(probepitch[in_range], yawcal2, yawcal1, period=360)

    #Dynamic calculations
    pitchbigger = (abs(probepitch)>abs(probeyaw))*1
    probemax = pitchbigger*probepitch + (pitchbigger -1)*probeyaw
    dyn[in_range] = np.interp(probemax[in_range], yawcal1, dyncal)

    #Find and resolve velocities
    U = np.sqrt(2 *-dyn * np.mean(BarnData.iloc[:, 0:4], axis = 1)/rho)
    U = np.real(U)


    #Resolve velocities in x, y and z directions
    Ux = U*np.cos(pitch*np.pi/180)*np.cos(yaw*np.pi/180)
    Uy = U*np.cos(pitch*np.pi/180)*np.sin(yaw*np.pi/180)
    Uz = U*np.sin(pitch*np.pi/180)
    
    

    return(yaw,pitch,U,Ux,Uy,Uz)

def findlongesttide(arr, min_length):
    # Find non-NaN values in the array
    not_nans = ~np.isnan(arr)
    
    # Compute differences between adjacent non-NaN values
    diffs = np.diff(not_nans.astype(int))
    
    # Find the start and end indices of each sequence without NaNs
    start_idxs = np.where(diffs == 1)[0] + 1
    end_idxs = np.where(diffs == -1)[0]
    
    # Handle the case where the first or last value is not NaN
    if not_nans[0]:
        start_idxs = np.concatenate([[0], start_idxs])
    if not_nans[-1]:
        end_idxs = np.concatenate([end_idxs, [len(arr)-1]])
    
    # Compute the length of each sequence without NaNs
    lengths = end_idxs - start_idxs + 1
    # Find the indices of sequences without NaNs that meet the minimum length
    long_idxs = np.where(lengths >= min_length)[0]

    longest_idx = np.argmax(lengths)
    
    if len(long_idxs) == 0:
        # No sequences meet the minimum length, return None
        return ([0, len(arr)], 0, len(arr))
    else:
        return ([(start_idxs[i], end_idxs[i]) for i in long_idxs], start_idxs[longest_idx], end_idxs[longest_idx])

@st.cache_data
def findintensity(U, Ux, Uy, Uz):

    #Find Turbulence Intensity in all dimensions
    intensity = (np.std(U))/(np.mean(U))*100
    intensityx = (np.std(Ux))/(np.mean(U))*100
    intensityy = (np.std(Uy))/(np.mean(U))*100
    intensityz = (np.std(Uz))/(np.mean(U))*100
        
    return(intensity, intensityx, intensityy, intensityz)

@st.cache_data
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

def callback():
    #Button was clicked
    st.session_state.button_clicked = not st.session_state.button_clicked
    st.session_state.low = low
    st.session_state.high = high

@st.cache_data
def datatocsv(data):
    return pd.DataFrame(index = False, sep = ',').encode('utf-8')
     

# Set the page layout and title
st.set_page_config(page_title='Barnacle Dashboard', page_icon='🌊', layout='wide', initial_sidebar_state= 'collapsed')
st.title('Barnacle Dashboard 📈')
# Define the chart layout
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


tab1, tab2, tab3 = st.tabs(['Upload Files', 'Data Presentation', 'Help'])

low, high = [], []
traces_U, traces_ux, traces_uy, traces_uz, traces_pitch, traces_yaw = [], [], [], [], [], []

with tab2:
    placeholderuploaddata = st.empty()

if "button_clicked" not in st.session_state:
        st.session_state.button_clicked = False

if "tide" not in st.session_state:
    st.session_state.tide = False

if "alltides" not in st.session_state:
    st.session_state.alltides = False

if "legend" not in st.session_state:
    st.session_state.legend = False

fig_pitch = []
fig_yaw = []

uploadedBarnFiles = []
numberoffiles = 0
i = 0
length_scales = []
with tab1:        
        
        col11, col12, col13 = st.columns(3)
        with col12:
            #st.write("Upload a file with the Barnacle Data and Zero Flow Data.")
            
            ZerosFile = st.file_uploader("Upload File(s) with Zeros Data", accept_multiple_files=  True)       
            
        
            if not ZerosFile:
        
                st.subheader('No Zero Flow file uploaded')
                st.session_state.button_clicked =  False
                with tab2:
                    with placeholderuploaddata.container():
                        st.subheader("Upload some Data to get started")

            if ZerosFile:
                for ZeroFile in ZerosFile:
                    ZerosData = pd.read_table(ZeroFile, sep = ',')
                    container = st.container()
                    st.subheader(f"Upload Barnacle Data associated with the {ZeroFile.name} Zero File")
                    BarnFiles = st.file_uploader("Upload File(s) with Barnacle Data", accept_multiple_files= True, key = f"BarnFilesUpload{ZeroFile}")
                    
                    if not BarnFiles:
                        st.subheader('No Barnacle Data uploaded')
                        with tab2:
                            with placeholderuploaddata.container():
                                st.subheader("Upload some Data to get started")
                    
                    with container:
                        if st.session_state.button_clicked == False:
                            low, high = st.slider(f"Select the region of the Zero flow file with genuine slackwater for {ZeroFile.name}:", 0, len(ZerosData), (0, len(ZerosData)), key = ZeroFile)
                            st.button('Confirm the Slack Water Time are correct', on_click=callback, use_container_width= True, key = f"confirm{ZeroFile}", type = "primary")
                        else:
                            st.button('Set Zero Flow region', on_click= callback,  use_container_width= True, key = f"goback{ZeroFile}")
                            
                            
                    if BarnFiles and ZerosFile:
                        
                        

                        if st.session_state.button_clicked:  
                            
                                for BarnFile in BarnFiles:
                                        
                                        uploadedBarnFiles.append(BarnFile)
                                        numberoffiles = numberoffiles + 1
                                        #Importing Data from text files
                                        
                                        BarnData = pd.read_table(BarnFile, sep = ',')
                                        yaw,pitch,U,Ux,Uy,Uz = DataProcessing(BarnData, ZerosData, st.session_state.low, st.session_state.high)
                                        tide_ranges, maxstart, maxend = findlongesttide(U, 5000)

                                        stacked_data = np.column_stack((yaw, pitch, U, Ux, Uy, Uz))

                                        # Create a DataFrame with appropriate column names
                                        columns = ['Yaw', 'Pitch', 'U', 'Ux', 'Uy', 'Uz']
                                        data = pd.DataFrame(stacked_data, columns=columns)

                                        csv = datatocsv(data)

                                        st.download_button(label = "Download Processed Data", data = csv, file_name = "data", mime = 'text/csv', use_container_width= True)
                                        intensities, intensitiesx, intensitiesy, intensitiesz = [], [], [], []
                                        if len(tide_ranges) > 3:
                                            for start,end in tide_ranges:
                                                eachintensity, eachintensityx, eachintensityy, eachintensityz = findintensity(U[start:end], Ux[start:end], Uy[start:end], Uz[start:end])
                                                intensities.append(eachintensity)
                                                intensitiesx.append(eachintensityx)
                                                intensitiesy.append(eachintensityy)
                                                intensitiesz.append(eachintensityz)
                                                intensity = round(np.mean(intensities),2)
                                                intensityx = round(np.mean(intensitiesx),2)
                                                intensityy = round(np.mean(intensitiesy),2)
                                                intensityz = round(np.mean(intensitiesz),2)

                                                alllength_scale = findlengthscale(U[start:end])
                                                length_scale = np.mean(alllength_scale)                                       
                                                length_scale = round(length_scale,2)
                                        else:
                                                start,end = tide_ranges
                                                eachintensity, eachintensityx, eachintensityy, eachintensityz = findintensity(U[start:end], Ux[start:end], Uy[start:end], Uz[start:end])
                                                intensities.append(eachintensity)
                                                intensitiesx.append(eachintensityx)
                                                intensitiesy.append(eachintensityy)
                                                intensitiesz.append(eachintensityz)
                                                intensity = round(np.mean(intensities),2)
                                                intensityx = round(np.mean(intensitiesx),2)
                                                intensityy = round(np.mean(intensitiesy),2)
                                                intensityz = round(np.mean(intensitiesz),2)

                                                alllength_scale = findlengthscale(U[start:end])
                                                length_scale = np.mean(alllength_scale)                                       
                                                length_scale = round(length_scale,2)
                                                length_scales.append(length_scale)
                                        times = np.arange(len(Ux)) / fs
                                        
                                        
                                        if st.session_state.tide:
                                            st.session_state.alltides = False
                                            #Create Velocity Traces for longest tides
                                            traces_U.append(go.Scatter(x = times, y = U[maxstart:maxend], name = f"{BarnFile.name} U - Intensity = {intensity}%"))
                                            traces_ux.append(go.Scatter(x = times, y = Ux[maxstart:maxend], name = f"{BarnFile.name} Ux - Intensity = {intensityx}%"))
                                            traces_uy.append(go.Scatter(x = times, y = Uy[maxstart:maxend], name = f"{BarnFile.name} Uy - Intensity = {intensityy}%"))
                                            traces_uz.append(go.Scatter(x = times, y = Uz[maxstart:maxend], name = f"{BarnFile.name} Uz - Intensity = {intensityz}%"))
                                            traces_pitch.append(go.Scatter(x = times, y=pitch[maxstart:maxend], name=f'{BarnFile.name} Pitch', mode='lines'))
                                            traces_yaw.append(go.Scatter(x = times, y=yaw[maxstart:maxend], name=f'{BarnFile.name}Yaw', mode='lines'))

                                        else:
                                            if st.session_state.alltides: 
                                                st.session_state.tide = False
                                                st.session_state.legend = False
                                                if len(tide_ranges) > 3:
                                                    for start,end in tide_ranges:
                                                        #Create Velocity Traces for longest tides

                                                            traces_U.append(go.Scatter(x = times, y = U[start:end], opacity = 0.75,  name = f"{BarnFile.name} U - Intensity = {intensity}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_ux.append(go.Scatter(x = times, y = Ux[start:end], opacity = 0.75, name = f"{BarnFile.name} Ux - Intensity = {intensityx}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_uy.append(go.Scatter(x = times, y = Uy[start:end], opacity = 0.75, name = f"{BarnFile.name} Uy - Intensity = {intensityy}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_uz.append(go.Scatter(x = times, y = Uz[start:end], opacity = 0.75, name = f"{BarnFile.name} Uz - Intensity = {intensityz}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_pitch.append(go.Scatter(x = times, y=pitch[start:end], opacity = 0.75, name=f'{BarnFile.name} Pitch, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}', mode='lines'))
                                                            traces_yaw.append(go.Scatter(x = times, y=yaw[start:end], opacity = 0.75, name=f'{BarnFile.name}Yaw, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}', mode='lines'))

                                                    else:
                                                            start,end = tide_ranges
                                                            traces_U.append(go.Scatter(x = times, y = U[start:end], opacity = 0.75,  name = f"{BarnFile.name} U - Intensity = {intensity}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_ux.append(go.Scatter(x = times, y = Ux[start:end], opacity = 0.75, name = f"{BarnFile.name} Ux - Intensity = {intensityx}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_uy.append(go.Scatter(x = times, y = Uy[start:end], opacity = 0.75, name = f"{BarnFile.name} Uy - Intensity = {intensityy}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_uz.append(go.Scatter(x = times, y = Uz[start:end], opacity = 0.75, name = f"{BarnFile.name} Uz - Intensity = {intensityz}%, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}"))
                                                            traces_pitch.append(go.Scatter(x = times, y=pitch[start:end], opacity = 0.75, name=f'{BarnFile.name} Pitch, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}', mode='lines'))
                                                            traces_yaw.append(go.Scatter(x = times, y=yaw[start:end], opacity = 0.75, name=f'{BarnFile.name}Yaw, time = {datetime.timedelta(seconds = int(start))}:{datetime.timedelta(seconds = int(end))}', mode='lines'))

                                            else:
                                                #Create Velocity Traces for full data set
                                                traces_U.append(go.Scatter(x = times, y = U, name = f"{BarnFile.name} U - Intensity = {intensity}%"))
                                                traces_ux.append(go.Scatter(x = times, y = Ux, name = f"{BarnFile.name} Ux - Intensity = {intensityx}%"))
                                                traces_uy.append(go.Scatter(x = times, y = Uy, name = f"{BarnFile.name} Uy - Intensity = {intensityy}%"))
                                                traces_uz.append(go.Scatter(x = times, y = Uz, name = f"{BarnFile.name} Uz - Intensity = {intensityz}%"))
                                                traces_pitch.append(go.Scatter(x = times, y=pitch, name=f'{BarnFile.name} Pitch', mode='lines'))
                                                traces_yaw.append(go.Scatter(x = times, y=yaw, name=f'{BarnFile.name}Yaw', mode='lines'))
                                                            
                                        #Create plots for yaw and pitch
                                        
                                        fig_pitch = go.Figure(data = traces_pitch, layout=chart_layoutpitch)
                                                                
                                        
                                        fig_yaw = go.Figure(data = traces_yaw, layout=chart_layoutyaw)
                                        
                                                                             

                    st.subheader("Navigate to Data Presentation tab")
                
            with tab2:
                if numberoffiles == len(uploadedBarnFiles):
                    i = 0
                    if 'length_scales' in locals():
                        for BarnFile in uploadedBarnFiles:
                            st.subheader(f"Estimated Length Scale for {BarnFile.name} data set is {length_scales[i]}m")
                            i = i+1
                        
                        # Display the charts
                        col1, col2 = st.columns(2, gap = 'medium')

                        with col1: 
                            toplot = st.multiselect('Select Velocity Component(s) to display', options = ['Velocity Magnitude', 'Ux', 'Uy', 'Uz'])
                            
                            traces_vel = []
                            for option in toplot:
                                if option == 'Velocity Magnitude' :
                                    traces_vel.extend(traces_U)
                                elif option == 'Ux':
                                    traces_vel.extend(traces_ux)
                                elif option == 'Uy':
                                    traces_vel.extend(traces_uy)
                                elif option == 'Uz':
                                    traces_vel.extend(traces_uz)
                            
                            fig_vel = go.FigureWidget(data=traces_vel, layout=chart_layout)
                            fig_vel.update_layout(showlegend = True)
                            st.plotly_chart(fig_vel, use_container_width=True)

                        with col2:
                            one, two = st.columns(2)
                            with one:
                                st.session_state.tide = st.checkbox(label = "Plot only the longest tide")
                            with two:
                                st.session_state.alltides = st.checkbox(label = "Plot all tides")
                            if fig_pitch:
                                st.plotly_chart(fig_pitch, use_container_width= True)
                            if fig_yaw:
                                st.plotly_chart(fig_yaw,  use_container_width= True)
                                    

with tab3:
    col1, col2, col3 = st.columns(3, gap = 'medium')
    with col2:
        st.header('Help')
        st.write('---')

        st.subheader('For pre-existing Barnacle Data')
        st.write('Use the Upload Data tab to upload Zero-Flow Data, multiple files can be uploaded. Upload the Barnacle Data and set the region of genuine slackwater for each Zero-Flow file and confirm. A button to download the processed data as a csv is at the bottom of the page.')
        st.write('Move to the Data Presentation tab to view the files that have been uploaded, the length scale for each data file is shown at the top of the page. Multiple velocity components can be shown on the graph and intensity can be seen in the legend.')
        st.write('Plots can be saved and manipulated using the tools in the top right of each plot')
        
        st.subheader('For Live Data')
        st.write('Select the Live Data tab using the sidebar.')
        st.write('Upload a Zeros flow file and set the region of slack water, confirm and press the data button to simulate live Barnacle Data.')

        st.subheader('Programming')
        st.write('Use the Programming tab to set the sampling rate of the Barnacle, and the expected depth and temperature of operation.')


