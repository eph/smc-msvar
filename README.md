Companion Code for 
================================================================
by MARK BOGNANNI and ED HERBST 


Requirements
------------
You need a 64-bit linux installation. 


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
