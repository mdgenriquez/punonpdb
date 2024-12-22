import streamlit as st
import pandas as pd
import requests
import io  # Importar el módulo io

# Título de la aplicación
st.title("Visualización de Moléculas (Sin RDKit)")

# URL del archivo
github_url = "https://raw.githubusercontent.com/mdgenriquez/punonpdb/main/MDG202_cl2.csv"

try:
    response = requests.get(github_url)
    response.raise_for_status()
    df = pd.read_csv(io.StringIO(response.text))  # Usar io.StringIO para leer el contenido

    if 'SMILES' not
