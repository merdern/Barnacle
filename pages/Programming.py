import streamlit as st

#Define page layout and title
st.set_page_config(page_title='Barnacle Dashboard', page_icon='🌊', layout='centered', initial_sidebar_state= 'expanded')
st.title('Barnacle Dashboard 📈')

tab1, tab2 = st.tabs(['Program Barnacle', 'Generate Zero Flow File'])

with tab1:

  st.header('Program the Barnacle')
  
  freq = st.number_input('Set Frequency of Data Logger (Hz)', value = 16)

  depth = st.number_input('Set Expected Depth of Barnacle (m)', value = 100)

  temp = st.number_input('Set Expected Temperature of Barnacle (°C)', value = 12)

  st.button("Program Barnacle", use_container_width=True)

with tab2:
  
  st.header('Generate a Zero Flow File')
  
  st.write('Place the Barnacle in slackwater and connect to the PC using a USB.')
  
  st.write('When the Barnacle is found by the software a generate zero-file button will appear. Press this to create a zero-file for the Barnacle. A download button will appear to save this.')
  
  st.button('Generate zero-flow file', use_container_width = True)
  
  st.button('Save zero-flow file', user_container_width = True)
