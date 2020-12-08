# set base image
FROM tensorflow/tensorflow:2.2.0

#RUN apt-get update && apt-get install --no-recommands -y \
#   package-to-install \
# && rm -rf /var/lib/apt/lists/*

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
