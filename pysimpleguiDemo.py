import PySimpleGUIWeb as sg
import datetime
import os
import numpy as np
import astropy.io.fits as fits
import glob
"""
  Demonstration of running PySimpleGUI code in repl.it!
  
  This demo program shows all of the PySimpleGUI Elements that are available for use.  New ones are being added daily.
  
  Now you can run your PySimpleGUI code in these ways:
  1. tkinter
  2. Qt (pyside2)
  3. WxPython
  4. Web Browser (Remi)
  5. repl.it (Remi)

  You can use repl.it to develop, test and share your code.
  If you want to run your GUI on tkinter, then change the import statement to "import PySimpleGUI".  To run it on WxPython, change it to "import PySimpleGUIWx".

  repl.it opens up an entirely new way of demonstrating problems, solutions, bugs, etc, in a way that doesn't require anything but a web browser.  No need to install a GUI package like tkinter.  No need to install Python for that matter.  Just open the repl link and have fun.

"""
print('Starting up...')

sg.ChangeLookAndFeel('LightGreen')      # set the overall color scheme
files = glob.glob('BAL_Plots/spectra_overplot_Plaw_*.png')
master = np.genfromtxt('Master_full_sample_redshift.txt',names=['name','ra','dec','zvi','tag'],dtype=('|S25',float,float,float,'|S10'))
pnames = [str(os.path.basename(sname).split('.')[0].split('_')[-1]).encode() for sname in files]
masternames=[sname for sname in master['name'] if sname in pnames]
masterra =[] ; masterdec=[];masterzvi = []
for name in masternames:
    yy=np.where( master['name']==name)[0][0]
    masterra.append(master['ra'][yy])
    masterdec.append(master['dec'][yy])
    masterzvi.append(master['zvi'][yy])

masternames= [sname.decode() for sname in masternames]
#masterra = [sra for sra in master['ra']]
#masterdec = [sdec for sdec in master['dec']]
#masterzvi = [szvi for szvi in master['zvi']]
mag = fits.open('FullSample3028Objs_query_for_SDSSMags.fits')[1].data

masterrmag = np.zeros(len(masternames))
for i in range(len(masternames)):
    xx=np.where(mag['Name'] == masternames[i])[0]
    if len(xx) >0:
        masterrmag[i] = mag['modelMag_r'][xx[0]]

masterrmag = ['{0:3.3f}'.format(srmag) for srmag in masterrmag]

# The GUI layout
layout =  [
        [sg.Text('',size=(50,1)),sg.Text('',size=(10,1)),sg.Text('SDSS BAL Factory: Visual Inspection', size=(140,2), font=('Comic sans ms', 40), text_color='red')],
        [sg.Text('',size=(50,1)),sg.Image(filename='/Users/vzm83/SDSSIV_BALQSO/pngtree-creative-simple-dividing-line.png',visible=True)],
        [sg.Text('',size=(30,1)),sg.Button('Previous Quasar', font=('Comic sans ms', 20)),sg.Text('', size=(5,1), key='_TEXT_008'), sg.Button('Next Quasar', font=('Comic sans ms', 20)),sg.Text('', size=(25,1), key='_TEXT_009'),sg.Button('MARK for future', font=('Comic sans ms', 20)),sg.Text('', size=(25,1)),sg.Button('WRITE to Database', font=('Comic sans ms', 20)),sg.Text('', size=(25,1)),sg.Button('Reset', font=('Comic sans ms', 20)),sg.Text('', size=(20,1)),sg.Text('', size=(5,1)),sg.Button('Exit', font=('Comic sans ms', 20))],
        [sg.Text('',size=(30,1))],
        [sg.Text('',size=(30,1)),sg.Text('NAME: ', size=(10,1)),sg.Text('', size=(12,1), key='_NAME_'),sg.Text('RA: ', size=(10,1)),sg.Text('', size=(8,1), key='_RA_'),sg.Text('DEC: ', size=(10,1)),sg.Text('', size=(8,1), key='_DEC_'),sg.Text('r-mag: ', size=(8,1)),sg.Text('', size=(8,1), key='_RMAG_'),sg.Text('Z_VI: ', size=(20,1)),sg.Text('', size=(6,1), key='_ZVI_')],
            [sg.Text('',size=(10,1)),sg.Text(''),sg.Image(filename='/Users/vzm83/SDSSIV_BALQSO/BAL_Plots/spectra_overplot_Plaw_Target0001.png',size=(1800,900), visible=True,key='_FILE_')],#sg.Spin(values=('Target0001','Target0002','Target0003','Target0004','Target0005'),size=(150,150),enable_events=True, key='_LISTnames_')],
            #[sg.Text('', size=(30,1), key='_TEXT_000')],
            # [sg.MultilineOutput('Multiline Output', size=(40,8),  key='_MULTIOUT_', font='Courier 12')],
            #[sg.Output(font='Courier 11', size=(60,8))],
            [sg.Text('',size=(20,1)),sg.Text('BAL Class', size=(40,1), font=('Comic sans ms', 20), text_color='Blue')],
            [sg.Text('',size=(20,1)),sg.Checkbox('HiBAL', enable_events=True, default=True, key='_hibal_'), sg.Checkbox('LoBAL', enable_events=True, key='_lobal_'),sg.Checkbox('FeLoBAL', enable_events=True, key='_felobal_'),sg.Checkbox('NTroughs', enable_events=True, key='_manynarrowsys_'),sg.Checkbox('J0300Analog', enable_events=True, key='_j0300analog_'),sg.Checkbox('Reddened', enable_events=True, key='_reddened_'),sg.Checkbox('IronEmitter', enable_events=True, key='_ironemitter_'),sg.Checkbox('Redshifted', enable_events=True,size=(20,1), key='_redshifted_'),sg.Checkbox('Acceleration', enable_events=False,size=(20,1), key='_accn_'),sg.Checkbox('Stable', enable_events=False,size=(20,1), key='_stable_'),sg.Checkbox('Emergence', enable_events=False,size=(20,1), key='_emergence_'),sg.Checkbox('Disappearance', enable_events=False,size=(20,1), key='_disappear_'),sg.Text('#troughs'),sg.Spin(values=('1', '2', '3','4','5','6'),initial_value='1', size=(10,1), enable_events=True, key='_nTrough_')],
            [sg.Text('',size=(20,1)),sg.Text('Absorption lines list', size=(30,1), font=('Comic sans ms', 20), text_color='Blue')],
            [sg.Text('',size=(20,1)),sg.Checkbox('C IV', enable_events=True, default=True, key='_CIV_'), sg.Checkbox('Si IV', default=True, enable_events=True, key='_SiiV_'),sg.Checkbox('Al III', enable_events=True, key='_AlIII_'),sg.Checkbox('Mg II', enable_events=True, key='_MgII_'),sg.Checkbox('Fe II', enable_events=True, key='_FeII_'),sg.Checkbox('Fe III', enable_events=True, key='_FeIII_'),sg.Checkbox('Al II', enable_events=True, key='_AlII_'),sg.Checkbox('Lo2HiBAL', enable_events=False,size=(20,1), key='_lo2hibal_'),sg.Checkbox('Hi2LoBAL', enable_events=False,size=(20,1), key='_hi2lobal_'),sg.Checkbox('Lo2FeBAL', enable_events=False,size=(20,1), key='_lo2febal_'),sg.Checkbox('Fe2LoBAL', enable_events=False,size=(20,1), key='_fe2lobal_'),sg.Checkbox('Co-ord Var', enable_events=False,size=(20,1), key='_coordvar_'),sg.Checkbox('Xtreme Spec Var', enable_events=False,size=(20,1), key='_xtreme_')],            #[sg.Combo(values=['Combo 1', 'Combo 2', 'Combo 3'], default_value='Combo 2', key='_COMBO_',enable_events=True, readonly=False, tooltip='Combo box', disabled=False, size=(12,1))],
            [sg.Text('', size=(20,1), key='_TEXT_001')],
            [sg.Text('',size=(20,1)),sg.Text('Future Observations ?', size=(35,1), font=('Comic sans ms', 20), text_color='Orange'),sg.Text('Score', size=(18,1), font=('Comic sans ms', 20), text_color='Orange'),sg.Text('Notes about BAL variability', size=(40,1),font=('Comic sans ms', 20), text_color='Orange', key='_TEXT_005'),sg.Text('Notes about quasar in general',font=('Comic sans ms', 20), text_color='Orange', size=(30,1), key='_TEXT__004')],
            [sg.Text('',size=(20,1)),sg.Listbox(values=('Must get a follow up', 'May be interesting', 'Typical Source'), size=(20,4), enable_events=True, key='_LIST_'), sg.Text('', size=(10,1), key='_TEXT_002'),sg.Slider((1,10), default_value=1, key='_SLIDER_', visible=True, enable_events=True),sg.Multiline('', do_not_clear=True, size=(40,4), enable_events=True,key='_VAR_COMMENT_'),sg.Multiline('', do_not_clear=True, size=(40,4), enable_events=True,key='_QUASAR_COMMENT_'),sg.Text('http://skyserver.sdss.org/dr15/en/tools/explore/summary.aspx?ra=0.00591&dec=20.01226', size=(30,1), key='_SDSS_')]
            #[sg.Text('', size=(30,1), key='_TEXT_003')],
            #[sg.Text('Notes about quasar in general', size=(30,1), key='_TEXT__004')],
            #[sg.Input('', do_not_clear=True, enable_events=True, size=(70,1))],
           # [sg.Text('Notes about BAL variability', size=(30,1), key='_TEXT_005')],
            #[sg.Multiline('', do_not_clear=True, size=(10,4), enable_events=True)],
            #[sg.Text('', size=(30,1), key='_TEXT_006')],
            #[sg.Button('Write to database'), sg.Text('', size=(10,1), key='_TEXT_007'),sg.Button('Mark for further inspection')],
            #[sg.Spin(values=('Write','Wait'),initial_value='Wait', size=(5,2))],
            #[sg.Button('Previous Quasar'),sg.Text('', size=(5,1), key='_TEXT_008'), sg.Button('Next Quasar'),sg.Text('', size=(5,1), key='_TEXT_009'),sg.Button('Reset'),sg.Text('', size=(5,1)),sg.Button('Exit')]
          ]

# create the "Window"
window = sg.Window('Visual Inspection of SDSS BAL QSO spectra',
                  layout=layout,
                   default_element_size=(12,1),
                   font='Helvetica 18',
                   )

start_time = datetime.datetime.now()
counter =0
def readwindow(counter=0):
    #event, values = window.Read(timeout=10)
    print(masternames[counter],type(masternames[counter])) 
    window.Element('_NAME_').Update(masternames[counter])
    window.Element('_RA_').Update(str(masterra[counter]))
    window.Element('_DEC_').Update(str(masterdec[counter]))
    window.Element('_RMAG_').Update(str(masterrmag[counter]))
    window.Element('_ZVI_').Update(str(masterzvi[counter]))
    print('Check-I')
    window.Element('_FILE_').Update('/Users/vzm83/SDSSIV_BALQSO/BAL_Plots/spectra_overplot_Plaw_{}.png'.format(masternames[counter]))
    print('Check-II')
    window.Element('_SDSS_').Update('http://skyserver.sdss.org/dr15/en/tools/explore/summary.aspx?ra={}&dec={}'.format(masterra[counter],masterdec[counter]))
    print('Check-III')
    #window.Element('_QUASAR_COMMENT_').Update('')
    #window.Element('_VAR_COMMENT_').Update('')
    npzfile = 'VI_out/BALfactory_VI_comments_{}.npz'.format(masternames[counter])
    if os.path.isfile(npzfile):
        ndata= dict(np.load(npzfile))
        print(ndata)
        for keys,vals in ndata.items():
            print(keys,vals,type(keys),type(vals))
            if ((keys == '_QUASAR_COMMENT_') |(keys=='_VAR_COMMENT_')):
                #continue
                strrep = vals.tobytes().decode()
                window.Element(keys).Update(str(vals))
            elif (keys != '_LIST_'):
                print(keys,vals)
                window.Element(keys).Update(vals)
    else:
        resetwindow(counter)
        resetwindow(counter)
    print('Check-IV')
    return event,values

def resetwindow(counter=0):
    window.Element('_QUASAR_COMMENT_').Update('')
    window.Element('_VAR_COMMENT_').Update('')
    window.Element('_SLIDER_').Update(1)
    #window.Element('_LIST_').Update('Typical Source')
    window.Element('_hibal_').Update(True)
    window.Element('_lobal_').Update(False)
    window.Element('_felobal_').Update(False)
    window.Element('_manynarrowsys_').Update(False)
    window.Element('_j0300analog_').Update(False)
    window.Element('_reddened_').Update(False)
    window.Element('_ironemitter_').Update(False)
    window.Element('_redshifted_').Update(False)
    window.Element('_nTrough_').Update(1)
    window.Element('_stable_').Update(False)
    window.Element('_accn_').Update(False)
    window.Element('_emergence_').Update(False)
    window.Element('_disappear_').Update(False)
    window.Element('_lo2hibal_').Update(False)
    window.Element('_hi2lobal_').Update(False)
    window.Element('_lo2febal_').Update(False)
    window.Element('_fe2lobal_').Update(False)
    window.Element('_xtreme_').Update(False)
    window.Element('_coordvar_').Update(False)
    window.Element('_CIV_').Update(True)
    window.Element('_SiiV_').Update(True)
    window.Element('_AlIII_').Update(False)
    window.Element('_MgII_').Update(False)
    window.Element('_AlII_').Update(False)
    window.Element('_FeII_').Update(False)
    window.Element('_FeIII_').Update(False)

    #return event,values




#  The "Event loop" where all events are read and processed (button clicks, etc)
while True:
    event, values = window.Read(timeout=10)     # read with a timeout of 10 ms
    if counter == 0:
        window.Element('_NAME_').Update(masternames[counter])
        window.Element('_RA_').Update(str(masterra[counter]))
        window.Element('_DEC_').Update(str(masterdec[counter]))
        window.Element('_RMAG_').Update(str(masterrmag[counter]))
        window.Element('_ZVI_').Update(str(masterzvi[counter]))
        #window.Element('_FILE_').Update('/Users/vzm83/SDSSIV_BALQSO/BAL_Plots/spectra_overplot_Plaw_{}.png'.format(masternames[counter]))


    #event,values = readwindow(counter)
    if ((event == 'Previous Quasar') & (counter >0)):
        event,values = readwindow(counter-1)
        event,values = readwindow(counter-1)
        counter -=1
    if ((event == 'Next Quasar')):
        event,values = readwindow(counter+1)
        event,values = readwindow(counter+1)
        counter +=1
    if event =='WRITE to Database':
        print('INSIDE WRITE LOOP')
        outsvzfile = 'VI_out/BALfactory_VI_comments_{}.npz'.format(masternames[counter])
        print(outsvzfile)
        if os.path.isfile(outsvzfile):
            mvcmd = 'mv {} VI_out/dump/.'.format(outsvzfile)
            print('Testing the MV cmd',mvcmd)
            os.system(mvcmd)
            np.savez(outsvzfile,**values)
        else:
            np.savez(outsvzfile,**values)
    if event == 'Reset':
        resetwindow(counter)
        resetwindow(counter)
        
    if event == 'MARK for future':
        print('INSIDE WRITE LOOP')
        futurefile = open('SDSSBALfactory_Marked_for_future.txt','a')
        futuretxt = '{}\t{}\t{}\n'.format(masternames[counter],masterra[counter],masterdec[counter])
        print(futuretxt)
        futurefile.write(futuretxt)
        futurefile.close()
    if event != sg.TIMEOUT_KEY:                 # if got a real event, print the info
        print(event, values)
        print('--------------------------------------------------------------------------------------------------')

        # also output the information into a scrolling box in the window
        # window.Element('_MULTIOUT_').Update(str(event) + '\n' + str(values), append=True)
    # if the "Exit" button is clicked or window is closed then exit the event loop
    if event in (None, 'Exit'):
        break
    #if event == 'Next Quasar':
        #window.Element('_DATE_').Update(str(datetime.datetime.now()-start_time))

    # Output the "uptime" statistic to a text field in the window
    #window.Element('_DATE_').Update(str(datetime.datetime.now()-start_time))

# Exiting the program
window.Close()    # be sure and close the window before trying to exit the program
print('Completed shutdown')

