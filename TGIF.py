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
    - make mp4 function work better
    - make this file prettier :)
    - would be pretty handy to have lightcurve option next to gif like in tesstransient with the little window moving along the lc
    - make it do tessreduce if no trobj exists
    - make multiple gifs/vids for an obs_list
    """
    
    def __init__(self, trobj=None,flux=None, mjd=None, time_window = 5000, size=90,
                 detection_list = None, objid=None, ra=None, dec=None, filepath=None, filename=None):
        """
        trobj: TESSreduce object.
        mjd: time of the event of interest
        time_window: the time in seconds to search before and after the specified MJD
        size: size of the TESSreduce cutout (not currently used)
        filepath: file path where gif/mp4 is saves
        filename: name of the gif/mp4
        """


        if flux is not None:
            self.flux=flux
        else:
            self.flux = self.trobj.flux

            

        
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

        if trobj == None:
            print('-- no tessreduce object given! --')
        else:
            self.trobj = trobj
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

    
        
    def make_frames(self,frames, vmin=0, vmax=20, fig_size=(4, 4),save_gif=False,save_mp4=False, dpi=100,interval=100, cmap="viridis", xlims=None, ylims=None):
        """
        make sure ffmpeg is installed properly or else mp4s will fail!
        """
        fig, ax = plt.subplots(figsize=fig_size)
        im = ax.imshow(self.flux[frames[0]], vmin=vmin, vmax=vmax, origin='lower', cmap=cmap)
        title = ax.set_title(f"frame {frames[0]}")
        ax.axis('off')

        if xlims is not None:
            ax.set_xlim(xlims)
        if ylims is not None:
            ax.set_ylim(ylims)

        def update(ind):
            im.set_array(self.flux[ind])
            title.set_text(f"frame {ind}")
            return im, title

        ani = animation.FuncAnimation(fig, update, frames=frames, interval=interval, blit=True)

        plt.show()
        if save_gif:
            ani.save(f"{self.filepath}/{self.filename}.gif", writer="pillow", dpi=dpi)
            print(f'saved GIF to {self.filepath}/{self.filename}.gif')
            plt.close()

        if save_mp4:
            ani.save(f"{self.filepath}/{self.filename}.mp4", writer="ffmpeg", dpi=dpi)
            print(f'saved MP4 to {self.filepath}/{self.filename}.mp4')
            plt.close()



