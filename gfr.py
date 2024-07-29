#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import gfr_read as gr
import gfr_rtf_reader as grr
import glob
import functools

def gfr():
    """
    Goes through a DICOMDIR and checks the year, if the year is before 2020
    and we find a PDF it means we have an older version of data storage.
    if it is 2020 or later it means we have the newer version and looks for
    Image storage to get the data from a examination.
    """
    dataset, path = gr.read_dcm()

    dirloc = path.split('DICOMDIR')[0]



# PDF: 1.2.840.10008.5.1.4.1.1.104.1
# Image 1.2.840.10008.5.1.4.1.1.7

    for series in dataset[0x00041220]:
        if (0x0008002a) in series:
            year = int(series[0x0008002a].value[:4])
        if 'year' in locals():
            if (0x00041510) in series:
                if '1.2.840.10008.5.1.4.1.1.7' in series[0x00041510].value:
                    if year >= 2020:
                        p = functools.reduce(lambda a,b : a + '/' + b, series[0x00041500])
                        total_path = dirloc + p
                        f, _ = gr.read_dcm(total_path)
                        print('after 2020')
                        d = gr.get_gfrdata(f)
                    else:
                        p0 = series[0x00041500][:3]
                        p = functools.reduce(lambda a, b : a.replace('DICOM',
                            'REPORTS') + '/' + b, p0)
                        total_path = dirloc + p
                        report_path = glob.glob(total_path +'/*')[0]
                        report = grr.read_report(report_path)
                        print('before 2020')
                        grr.get_data(report)

                elif '1.2.840.10008.5.1.4.1.1.104.1' in series[0x00041510].value:
                    if year < 2020:
                        p0 = series[0x00041500][:3]
                        p = functools.reduce(lambda a, b : a.replace('DICOM',
                            'REPORTS') + '/' + b, p0)
                        total_path = dirloc + p
                        report_path = glob.glob(total_path +'/*')[0]
                        report = grr.read_report(report_path)
                        print('before 2020')
                        grr.get_data(report)

                        #                        p = functools.reduce(lambda a, b : a + '/' + b, series[0x00041500])
#                        total_path = dirloc + p
#                        f, _ = gr.read_dcm(total_path)
#                        report = grr.read_report(total_path)
#                        print('before 2020')
#                        gr.get_gfrdata(report)
                    
                