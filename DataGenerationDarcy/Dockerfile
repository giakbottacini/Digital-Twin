# Usa l'immagine base di FEniCS
FROM ghcr.io/scientificcomputing/fenics-gmsh:2024-05-30

# Installa pyDOE2 e altre librerie necessarie
RUN pip install pyDOE2 tqdm scikit-learn

# Setta la directory di lavoro dentro al contenitore
WORKDIR /root/shared

# Copia il codice del progetto nel contenitore (se necessario)
COPY . /root/shared

# Comando che verrà eseguito quando il contenitore si avvia
CMD ["bash"]
