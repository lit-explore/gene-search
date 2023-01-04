FROM tiangolo/uvicorn-gunicorn-fastapi:latest

MAINTAINER keith.hughitt@gmail.com

# install python deps
WORKDIR /app

ADD requirements.txt /app/requirements.txt
RUN pip install -r requirements.txt

# add main api code
ADD ./app /app
