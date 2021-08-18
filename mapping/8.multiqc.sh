#!/bin/bash

module load MultiQC/1.9-foss-2019b-Python-3.7.4
cd $BASEDIR
multiqc .
