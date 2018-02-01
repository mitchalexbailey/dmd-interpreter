FROM continuumio/miniconda

# We need gcc and build tools to complie uwsgi
RUN apt-get install -y linux-headers-amd64 build-essential libpcre3 libpcre3-dev

# We need pip to install uwsgi, django, etc. (sqlite is for the database)
RUN conda install -y pip
RUN conda install -y sqlite

# Setup some interpreter-specific values as environmental variables for documentation
ENV INTERPRETER_DIR /interpreter
ENV INTERPRETER_USER interpreter
ENV INTERPRETER_UID 1000
ENV INTERPRETER_PORT 8080

EXPOSE $INTERPRETER_PORT

# Install the app server and application dependencies
RUN pip install uwsgi

RUN echo this is a hack

# Add this source directory to the image
ADD . $INTERPRETER_DIR
RUN pip install -r $INTERPRETER_DIR/requirements.txt
RUN cd $INTERPRETER_DIR ; python manage.py collectstatic --no-input

RUN useradd --no-create-home --home-dir $INTERPRETER_DIR --uid $INTERPRETER_UID --user-group $INTERPRETER_USER

# Set ownership of application files to runtime user and switch to that account
RUN chown -R $INTERPRETER_USER $INTERPRETER_DIR

CMD uwsgi --chdir=$INTERPRETER_DIR --http=0.0.0.0:$INTERPRETER_PORT --ini=$INTERPRETER_DIR/docker-uwsgi.ini
