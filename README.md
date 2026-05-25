# P-n Junction and Diode Simulator
This program provides a visual representation of a P-n junction under different applied voltage, as well as the corresponding I-V curve. The program has an interactive slider component that enables the user to freely explore the influence of applied voltage on a p-n junction.

LINK to the webpage: https://p-n-junction-and-diode-simulator.onrender.com/

This webpage was deployed using a free service Render, so it needs a bit of time to load completely.

I believe this program will come in handy in Semiconductor classes in universities.

Screenshots:


<img width="1204" height="774" alt="Screenshot 2026-05-25 143301" src="https://github.com/user-attachments/assets/55e872ac-f6ca-4781-b525-cbec8d866bbd" />

<img width="1163" height="592" alt="Screenshot 2026-05-25 143321" src="https://github.com/user-attachments/assets/c9440ada-4a1d-4943-920a-27a21aedf1bf" />



## Repository conten![Uploading Screenshot 2026-05-25 143301.png…]()
ts:
- main.py : initial commit that is built upon Matplotlib's pyplot chart
- app.py : website app made with Dash library for easy access and share
  
## Constants
For this visualization, I used the doping concentarions of 10^17 and 10^16 for Na and Nd respectively. It provides an interesting contrast between the acceptor and donors. You can change it freely in the Constants at the very beginning of the program.

P-n junction operates at 300K (room temperature).

## Future Improvements
I plan to incorporate a couple of features in the nearest future:
 - Choice of material (Si, GaAs, Ge, etc.)
 - MOSFET simulator inntegration
 - Solar Cell mode

Let me know if you have any suggestions or feedback! Thank you
