#!/bin/sh
#############################################################
#                                                           #
# Contact:                                                  #
#   Maarten Pater (m.pater@amc.uva.nl / www.mirdesign.nl)   #
#   Bjorn Koel (b.f.koel@amc.uva.nl)                        #
#                                                           #
#############################################################
# This script is used to enable continues execution of      #
# shell scripts.                                            #
#############################################################

cd $DIR_Start
if [ "$1" = "Ask" ] 
then
    echo "Would you like to execute '${2}' to start the next phase in the pipeline?"
    read -p "(Press Y or N and enter.)" yn
    case $yn in
        [Yy]* ) sh ./$2;;
    esac
elif [ "$1" = "Continious" ] 
then
    sh ./$2
fi