Companion Code for A SEQUENTIAL MONTE CARLO APPROACH TO INFERENCE IN MULTIPLE-EQUATION MARKOV-SWITCHING MODELS
================================================================
by MARK BOGNANNI and ED HERBST [ed.herbst@gmail.com]


Requirements
------------
You need a 64-bit linux installation, GCC 5+, and MPI.  


Installation + Usage
--------------------
1. Download Anaconda

2. Install the and activate msvar environment 
   ```sh
   conda env create -f env.yaml
   source activate msvar
   ```

3. Run the estimation `create_models.py`
   ```sh
   python create_models.py [-m {1,2}] [-v {1,2,3,4,5}] {swz,rfb,rfb-hier}
   ```

Note
----
The code was refactored recently, so there may be bugs.  Please let me know. 
