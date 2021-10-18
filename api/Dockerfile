FROM tiangolo/uvicorn-gunicorn-fastapi:latest

MAINTAINER keith.hughitt@gmail.com

# install neovim
RUN apt-get update && apt-get install -y neovim python3-neovim 
RUN python3 -m pip install pynvim

# aliases
RUN echo 'alias l="ls -lah"' >> ~/.bashrc
RUN echo 'alias vi="nvim"' >> ~/.bashrc
RUN echo 'alias vim="nvim"' >> ~/.bashrc

# install python deps
WORKDIR /app

ADD requirements.txt /app/requirements.txt
RUN pip install -r requirements.txt

# add main api code
ADD ./app /app
ADD ./annot /app/annot
