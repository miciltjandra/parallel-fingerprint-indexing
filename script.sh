#!/bin/sh
LOCAL_PATH="/home/chills/parallel-fingerprint-indexing/"
SERVER_PATH="/home/m13514108/parallel-fingerprint-indexing"
SERVER_LOGIN="m13514108@167.205.32.236"

if [ "$1" = "ssh" ]
then
    ssh $SERVER_LOGIN
elif [ "$1" = "scps" ]
then
    if [ "$2" != "" ]
    then
        scp $2 $SERVER_LOGIN:$SERVER_PATH
    fi
elif [ "$1" = "scpl" ]
then
    if [ "$2" != "" ]
    then
        scp $SERVER_LOGIN:$SERVER_PATH/$2 $LOCAL_PATH
    fi
fi