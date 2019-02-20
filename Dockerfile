FROM google/cloud-sdk:slim

RUN apt-get update
RUN apt-get install -y unzip wget zlib1g-dev g++ make

COPY requirements.txt .

RUN pip install -r requirements.txt

# copying python script to container runtime working directory for google-v2 provider
COPY test_reg.py .

CMD [ "/bin/bash" ]
