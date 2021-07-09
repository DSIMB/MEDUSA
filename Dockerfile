# set base image
FROM tensorflow/tensorflow:2.3.0

LABEL program="MEDUSA"
LABEL description="A Deep Learning based protein flexibility prediction tool."
LABEL version="1.1"
LABEL maintainer="gabriel.cretin@u-paris.fr"

# set the working directory in the container
WORKDIR /MEDUSA

# copy the dependencies file to the working directory
COPY requirements.txt .

# install dependencies
RUN pip install -r requirements.txt

# copy the content of the local src directory to the working directory
COPY bin/launch_MEDUSA.sh .

# command to run on container start
ENTRYPOINT ["./launch_MEDUSA.sh"]
CMD ["--help"]
