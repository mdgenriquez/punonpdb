import streamlit as st
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Draw

# Título de la aplicación
st.title("Visualización de Moléculas en 2D")
st.markdown("Sube un archivo CSV que contenga moléculas en formato SMILES para visualizarlas en 2D.")

# Subir archivo CSV
uploaded_file = st.file_uploader("Cargar archivo CSV", type="csv")

if uploaded_file:
    # Leer el archivo CSV
    try:
        df = pd.read_csv(uploaded_file)
        if 'SMILES' not in df.columns:
            st.error("El archivo debe contener una columna llamada 'SMILES'.")
        else:
            st.success("Archivo cargado correctamente.")

            # Mostrar tabla con las moléculas
            st.subheader("Tabla de Moléculas")
            st.dataframe(df)

            # Generar imágenes de moléculas
            st.subheader("Visualización de Moléculas en 2D")

            for index, row in df.iterrows():
                smiles = row['SMILES']
                try:
                    mol = Chem.MolFromSmiles(smiles)
                    if mol:
                        img = Draw.MolToImage(mol, size=(300, 300))
                        st.image(img, caption=f"Molécula {index + 1}: {smiles}", use_column_width=False)
                    else:
                        st.warning(f"No se pudo generar la estructura para el SMILES: {smiles}")
                except Exception as e:
                    st.error(f"Error procesando la molécula en la fila {index + 1}: {e}")
    except Exception as e:
        st.error(f"Error al leer el archivo: {e}")
