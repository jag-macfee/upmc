# Unsteady Panel Method Computer

## Overview

This project implements panel method solutions for analyzing airfoil behavior in unsteady flow conditions. The code solves for circulation distributions on flat plate airfoils experiencing gust and turbulent flow using both steady and unsteady panel methods.

## Getting Started

### Cloning the Repository

If you don't have Git installed:

- **macOS**: Git is pre-installed. Verify by typing `git --version` in Terminal
- **Windows**: Download from [git-scm.com](https://git-scm.com/downloads)
- **Linux**: Install via package manager: `sudo apt-get install git`

To clone this repository:

1. Open Terminal (macOS/Linux) or Command Prompt (Windows)
2. Navigate to where you want to store the project:
    ```bash
    cd ~/Desktop
    ```
3. Clone the repository using HTTPS (no SSH key required):
    ```bash
    git clone https://github.com/jag-macfee/upmc.git
    ```
    Replace the URL with the actual repository URL. You can find this by clicking the green "Code" button on GitHub and copying the HTTPS URL.
4. Navigate into the project folder:
    ```bash
    cd upmc
    ```

**Note**: HTTPS cloning works without any setup. If you prefer SSH, you'll need to [set up SSH keys with GitHub](https://docs.github.com/en/authentication/connecting-to-github-with-ssh) first.

## Python Setup

### Installing Python

If you don't have Python installed:

- **macOS**: Python 3 is pre-installed. Verify by opening Terminal and typing `python3 --version`
- **Windows**: Download from [python.org](https://www.python.org/downloads/) and run the installer
- **Linux**: Install via package manager: `sudo apt-get install python3 python3-pip`

### Installing Required Packages

Open Terminal (macOS/Linux) or Command Prompt (Windows) and navigate to the project folder:

```bash
cd /path/to/upmc
```

Install all required packages using the requirements file:

```bash
pip3 install -r requirements.txt
```

Alternatively, if `pip3` doesn't work, try:

```bash
pip install -r requirements.txt
```

**Minimal Installation**: If you only want to run the basic panel method code, you can install just the essential packages:

```bash
pip3 install numpy matplotlib
```

## How to Run

### Steady Panel Method

1. Open Terminal and navigate to the steady panel method folder:

    ```bash
    cd steady-panel-method
    ```

2. Run the analysis:

    ```bash
    python3 main.py
    ```

3. A plot comparing numerical and theoretical solutions will appear.

### Unsteady Panel Method

1. Navigate to the unsteady panel method folder:

    ```bash
    cd unsteady-panel-method
    ```

2. Run the analysis:
    ```bash
    python3 main.py
    ```

### Troubleshooting

- If `python3` doesn't work, try `python` instead
- Make sure you're in the correct directory before running `main.py`
- Ensure all required packages are installed (see Python Setup section)
