# Eidolon

This software was written in Python 2.7.6 (32 bit) on a Windows 7 machine.<br/>
It relies heavily on NumPy, version 1.8.1, and on PILLOW (PIL), version 1.1.7.

It is a translation of software written in the Processing language 
(see https://processing.org/) by <a href="http://www.gestaltrevision.be/en/about-us/contact/all-contacts/45">Jan Koenderink</a>.

The original as well as a matlab version (by <a href="http://www.allpsych.uni-giessen.de/matteo/">Matteo Valsecchi</a>) can be found <a href="http://www.allpsych.uni-giessen.de/EidolonFactories/index.htm">here</a>.

A manual from the original software is available in the manual folder (<a href="https://github.com/gestaltrevision/Eidolon/blob/master/manual/Manual_1024x512_page.pdf">PDF</a>).

The Python version relies rather heavily on generators, since there are sometimes a lot of matrices involved at a given time, and I ran into memory problems doing it any other way.<br/>
There’s a small library called eidolon (surprise) which contains the routines. There are 4 files there: helpers.py, noise.py, picture.py and scalespaces.py. This is where the magic happens.<br/>
The file helpers.py contains a number of essential functions, picture.py is a class that contains all the basic stuff the picture does and the 2 others contain classes that provide the necessary generators. 

The file eidolon.py just sets some parameters and calls upon the examples contained in examples.py. You can run both files (there’s a test function in examples.py). The eidolon.py file more or less mimics Jan’s start file in his original Processing (JAVA) version.
