# DomeBuilder - A simple dome optimiser
=======================================

## Overview

This builds a hybrid Schwedler-Lamella dome structure and analyses it.

DomeBuilder has two dependencies from GitHub:

* For the frame analysis: [FJFranklin/PyNite](https://github.com/FJFranklin/PyNite) (a fork of Craig Brinck's [JWock82/PyNite](https://github.com/JWock82/PyNite))
* For the optimiser: [FJFranklin/BeesEtAl](https://github.com/FJFranklin/BeesEtAl)

If these are downloaded as source in the same directory as DomeBuilder, you'll need to set Python's path, e.g.:
```
export PYTHONPATH=$PWD/BeesEtAl:$PWD/PyNite
cd DomeBuilder
python3 test.py
```

# License

Unless otherwise stated, the code and examples here are
provided under the MIT License:

Copyright (c) 2019 Francis James Franklin

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
