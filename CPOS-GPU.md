# CPOS GPU cluster

## Configuration
- Hostname: `hpcf3.cpos.hku.hk` (access only within HKU network)
- 4 computing nodes (g03, g04, g05, g06), each with 4 x L40 GPU cards, 120 thread CPU and 1TB of ram Total 10TB of SSD storage
- Duration is from 31 Mar to 6 April (~ 1 week)
- Jobs submission via PBS scheduler
- Project user home folders (/home2) and project folder (/mnt/project) creation,  with 10TB shared quota assigned

## Login for job management
As this is cluster with the PBS job management system for many computing nodes, 
it is generally recommended to login to the login node and sub/view your jobs.

To login, use the command on your terminal:

```bash 
ssh USER@hpcf3.cpos.hku.hk
```

Also, on this node, you may create symbolic link (somehow, it doesn't support via computing node):

```bash
ln -s /mnt/project/ ~/hackathon
```

## Login to computing node
For this special exercise with a small number of users, we can directly login 
and run jobs on computing nodes.

- hpcf3-g03 (GPU 1 to 4): Shumin, Ruiyan, Mingze, Julia
- hpcf3-g04 (GPU 1 to 4): Tianjie, Lingyu, Lousia, Xianjie
- hpcf3-g05 (GPU 1 to 4): Jiamu, Kevin, Chrissy, Bingxue
- hpcf3-g06 (GPU 1 to 4): Fangxin, Yuanhua, Amy, Minghao 

There are multiple ways, including using ProxyJump, either directly go from 
terminal or configure your 

1. directly use terminal (choose your own username and computing node above):

```bash
ssh USER@hpcf3-g03 -J USER@hpcf3.cpos.hku.hk
```

2. configure your `.ssh/config` file on your laptop/desktop:

This is much easier for use VS code or Jupyter lab on web browser:
```
Host hpcf3-g03 hpcf3-g04 hpcf3-g05 hpcf3-g06
    User YOUR_USERNAME
    ProxyJump YOUR_USERNAME@hpcf3.cpos.hku.hk
    ServerAliveInterval 5
```

**Note**, change the `YOUR_USERNAME` to your username!!

Then you can simply use (choose your own username and computing node above):

```bash
ssh hpcf3-g03
```

## Interactive jobs on computing node

For interactive jobs, we assume that you have configured your `.ssh/config` file 
as mentioned above.

Then you can choose one of the following to use Jupyter lab:

1. VS Code:

Very easy, same as other standard computing server, except you have to type 
password twich, one for login-node, one for the specific computing node.

2. Web browser:

- step 1: login to the computing node
  ```
  # Optionally, create an environment with jupyter lab
  module load miniconda3/24.11.1_py312
  conda init
  conda create -n myenv python==3.11 jupyterlab
  conda activate myenv
  ```

- step 2: open jupyter lab
  ```
  # change to any port you like, as long as not used by others
  jupyter lab --port 7111 --no-browser --ip=0.0.0.0
  ```

- step 3: type the command by adding computing node:
  ```
  ssh -N -L PORT:COMP_NODE:PORT USER@hpcf3.cpos.hku.hk
  ```
  Please change $PORT, $COMP_NODE and $USER to your own setting. Then go to your
  local browser and open https://localhost:PORT

