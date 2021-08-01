<h1 align="center">
	<img alt="Circuit-breaker" title="COMTRADE-procesing" src="./src/public/vishnu-mohanan-yQpAaMsQzYE-unsplash.jpg" />
</h1>


# COMTRADE processing

![Badge](https://img.shields.io/badge/license-MIT-informational?style=for-the-bade)
![Badge](https://img.shields.io/badge/Octave-v6.2.0-0790C0?logo=Octave) ![Badge](https://img.shields.io/badge/status-in_development-red)

## ❗ DISCLAIMER ❗

This software is a school assignment and is provided as is. __It is not fit for use in real world applications__. Use at your own risk.

## About

This repository contains a power plant automated fault detection programming assignment for the discipline [PEA3411 - Introduction to Electric Systems Automation](https://uspdigital.usp.br/jupiterweb/obterDisciplina?sgldis=PEA3411&verdis=2).

This program consists of a script written in Octave that simulates an electrical automation system that receives an signals from a power line and detects electrical faults with digital signal processing. The input for this program is any COMTRADE file generated by an IED.

# Table of contents
- [COMTRADE processing](#comtrade-processing)
  - [❗ DISCLAIMER ❗](#-disclaimer-)
  - [About](#about)
- [Table of contents](#table-of-contents)
- [Status](#status)
- [Features](#features)
- [How to install and run this application](#how-to-install-and-run-this-application)
  - [Prerequisites](#prerequisites)
  - [Installing and running](#installing-and-running)
- [Techs utilized](#techs-utilized)
- [Author and Acknowledgements](#author-and-acknowledgements)
# Status
This project is a __work in progress__. 
# Features
- [ ] File reading
- [ ] Signal pre-processing
- [ ] Internal protection
- [ ] Analog low-pass filtering
- [ ] ADC digitalization
- [ ] Digital signal processing with FFT
# How to install and run this application
## Prerequisites
To run this code you will need the latest versions of [Octave](https://www.gnu.org/software/octave/index) and a COMTRADE file. I'm sorry that I'm not able to make any COMTRADE files available for use, but they are of sensitive nature and can only be used with authorization from their proprietors.
## Installing and running
```bash
# Clone this respository
$ git clone <https://github.com/Gabriel-f-r-bojikian/COMTRADE-processing>

# Move into the newly created folder
$ cd COMTRADE-processing
```

You have now succesfully cloned this application code to your machine. You should now add your COMTRADE files to the root of the directory that you just created. You just need to open your Octave application the corresponding folder and change the `filename` variable in the `COMTRADE_processor.m` file to match the COMTRADE file you are analyzing.

# Techs utilized
The following technologies were used for the completion of this application:
- ![Octave Badge](https://img.shields.io/badge/-Octave-0790C0?style=flat-square&logo=Octave&logoColor=ffffff&link=https://nodejs.org/en/)



# Author and Acknowledgements
This project was done from scratch by me, Gabriel Fernandes, with several lines of code borrowed from the lectures given by Eduardo Lorenzetti Pellini (elpellini@usp.br), Full Professor of the Energy and Automation department of the Polytechnic School of Engineering of the University of São Paulo.

The banner image credit goes to [Vishnu Mohanan](https://unsplash.com/@vishnumaiea).
The image source can be found in this [link](https://unsplash.com/photos/yQpAaMsQzYE).

My contact methods are:
[![Github Badge](https://img.shields.io/badge/-Gabriel-181717?style=flat-square&logo=github&logoColor=white&link=https://github.com/Gabriel-f-r-bojikian)](https://github.com/Gabriel-f-r-bojikian) [![Linkedin Badge](https://img.shields.io/badge/-Gabriel-blue?style=flat-square&logo=Linkedin&logoColor=white&link=https://www.linkedin.com/in/gabriel-fernandes-rosa-bojikian-688b84164/)](https://www.linkedin.com/in/gabriel-fernandes-rosa-bojikian-688b84164/) [![Gmail Badge](https://img.shields.io/badge/-gabriel.bojikian.dev@gmail.com-c14438?style=flat-square&logo=Gmail&logoColor=white&link=mailto:gabriel.bojikian.dev@gmail.com)](mailto:gabriel.bojikian.dev@gmail.com)