FROM node:latest

MAINTAINER keith.hughitt@gmail.com

RUN mkdir -p /home/node/app
WORKDIR /home/node/app

# aliases
RUN echo 'alias l="ls -lah"' >> ~/.bashrc

ADD app/package.json ./
ADD app/yarn.lock ./

RUN yarn install

COPY --chown=node:node app/ /home/node/app
