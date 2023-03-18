import math
import numpy as np
pi = np.pi

def simple_pitch_inches1(r_vector_inches, Pitch_inches): #r_vector and pitch in inches
    r_vector = [rr*0.0254 for rr in r_vector_inches]
    Pitch = Pitch_inches*0.0254
    return simple_pitch(r_vector, Pitch)

def simple_pitch_inches2(r_vector, Pitch_inches): #only pitch in inches
    Pitch = Pitch_inches*0.0254
    return simple_pitch(r_vector, Pitch)

def simple_pitch(r_vector, Pitch): #All values in meters
    Beta_dist = []
    for rr in r_vector:
        Beta_dist.append(math.degrees(math.atan(Pitch/(2*pi*rr))))
    return Beta_dist