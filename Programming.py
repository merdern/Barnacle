import streamlit as st

#Define page layout and title
st.set_page_config(page_title='Barnacle Dashboard', page_icon='🌊', layout='centered', initial_sidebar_state= 'expanded')
st.title('Barnacle Dashboard 📈')

freq = st.number_input('Set Frequency of Data Logger', value = 16)

depth = st.number_input('Set Expected Depth of Barnacle', value = 100)

st.button("Program Barnacle")

