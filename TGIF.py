import numpy as np
import matplotlib.pyplot as plt
import tessreduce as tr
import pandas as pd
from astropy.wcs import WCS as astropyWCS
from astropy.coordinates import SkyCoord
import matplotlib.animation as animation
import os



def nice_plots(change=True):
    plt.rc('font', family='serif')
    fig_width_pt = 244.0
    text_width_pt = 508.0

    inches_per_pt = 1.0/72.27
    golden_mean = (np.sqrt(5)-1.0)/2.0  
    fig_width = fig_width_pt*inches_per_pt*1.5

    fig_width_full = text_width_pt*inches_per_pt
    fig_height =fig_width*golden_mean    
    fig_size = [fig_width,fig_height] 

def to_seconds(val):
    return val*24*60*60

def to_days(val):
    return val/(24*60*60)

def mag(flux, zeropoint=20.44):
    return -2.5*np.log10(flux) + zeropoint

class tgif:
    """
    TGIF = T(ESS) GIF (^:
    Uses a tessreduce object to create a small gif or mp4 of your region of interest!
    to-do:
    - implement working on ozstar
    - make mp4 function work again
    - make this file prettier :)
    - would be pretty handy to have lightcurve option next to gif like in tesstransient with the little window moving along the lc
    - make it do tessreduce if no trobj exists
    - make multiple gifs/vids for an obs_list
    """
    def __init__(self, trobj=None, mjd=None, time_window = 5000, size=90,
                 detection_list = None, objid=None, ra=None, dec=None, filepath=None, filename=None):
        """
        trobj: TESSreduce object.
        mjd: time of the event of interest
        time_window: the time in seconds to search before and after the specified MJD
        size: size of the TESSreduce cutout (not currently used)
        filepath: file path where gif/mp4 is saves
        filename: name of the gif/mp4
        """
        if trobj == None:
            print('-- no tessreduce object given! --')
        else:
            self.trobj = trobj

        self.mjd = mjd
        self.time_window = time_window # in seconds
        self.detection_list = None
        self.objid = None
        self.event_info = None
        self.size = size
        self.ra = ra
        self.dec = dec
        self.detection_list = detection_list
        self.objid = objid
        self.detection_list = None
        
        self.detected_inds = np.where(((self.trobj.lc[0] >= (self.mjd - to_days(self.time_window)))) & 
                                      (self.trobj.lc[0] <=  (self.mjd + to_days(self.time_window))))[0]
        
        self.filename = filename
        self.filepath = filepath
        if self.filepath is None:
            self.filepath = os.getcwd()


    def get_event_info(self):
        """
        Gets event info (ra, dec, mjd, objid) from TESSELLATE detection csvs.
        """
        self.event_info = self.detection_list[self.detection_list['objid']==self.objid]
        
        
    def make_gif(self, vmin=0, vmax=20, fig_size=(4, 4), side_frames=5, save=False):
        start_idx = self.detected_inds[0] - side_frames      
        n_frames = len(self.detected_inds) + side_frames        
        

        fig, ax = plt.subplots(figsize=fig_size)

        im = ax.imshow(self.trobj.flux[start_idx], vmin=vmin, vmax=vmax, origin='lower', cmap='viridis')
        title = ax.set_title(f"frame {start_idx}")
        ax.axis('off')

        def update(frame):
            idx = start_idx + frame
            im.set_array(self.trobj.flux[idx])
            title.set_text(f"frame {idx}")
            return [im, title]

        ani = animation.FuncAnimation(fig, update, frames=n_frames, interval=300, blit=True)
        if save:
            file_name= self.filename
            ani.save(f"{self.filepath}/{file_name}.gif", writer="pillow", dpi=100)
            print(f'saved to {self.filepath}/{file_name}.gif')
            plt.close()

