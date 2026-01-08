FROM python:3.12.3

#installing ncbi to run the blastn
RUN apt-get update && \
    apt-get install -y ncbi-blast+ && \
    apt-get clean


#Also need to install pgmauv
RUN pip install pygenomeviz
RUN apt install progressivemauve

# Install mummer for nucmer
RUN apt-get update && \
    apt-get install -y mummer && \
    apt-get clean

# setting up the work directory
WORKDIR /app

#copying requirements
COPY requirements.txt /app/requirements.txt
RUN pip install --no-cache-dir -r /app/requirements.txt

#copying the origenomi cli
COPY README.md /app/
COPY pyproject.toml /app/
COPY src /app/src/

#copying the fastapi server for the same tool
COPY temp_anant/test/app /app/app

#copying the databases
COPY clustered_DB /app/clustered_DB/

RUN pip install -e .

# the pmauve is not maintained since 2015 so most of the modern container can't use it so i had to install it on the 
# the system and the system's pmauve is used here not the one in the docker
# ENV PMAUVE_PATH=/mauve/progressiveMauve


CMD ["uvicorn", "app.main:app", "--host", "0.0.0.0", "--port", "8000"]


