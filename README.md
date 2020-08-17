# Distributed and parallel image processing

Please see the [proposal](tp2.pdf), [report](tp2-Relatório.pdf) or [presentation](tp2-Apresentação.pdf) for all the details.

## HellfireOS Realtime Operating System

---
### The operating system structure is organized as follows:

- app/ - applications
- arch/ - architecture specific device drivers and kernel bootstrapping
- drivers/ - general purpose device drivers
- fs/ - filesystem drivers
- lib/ - standard system / application libraries
- net/ - lightweight network stack
- platform/ - image building scripts for different platforms (application, kernel, architecture)
- sys/ - kernel core
- usr/ - simulators and documentation
