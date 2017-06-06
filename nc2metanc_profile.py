#!/usr/bin/python2.7

#  Combines NetCDF files containing separate CTD profiles into
#  a single NetCDF file containing many profiles.

import sys,os
import numpy as np
import netCDF4
import glob
import subprocess
import pyratemp
import requests,bs4

#  THE URL OF THE MACHINE RUNNING THE ERDDAP SERVER
#thredds_url = "http://gcoos4.tamu.edu:8080"
thredds_url = "http://138.197.0.26:8080"

#  THE URL OF THE GOMRI DATA SITE
gomri_url = "https://data.gulfresearchinitiative.org/data/"

#  THE FULL PATHNAME OF THE FOLDER CONTAINING THE ERDDAP NETCDF METAFILES
erddap_meta_filedir = "/data/erddap/profile/"

#  TEMPORARY PREFIX FOR ERDDAP TITLES (IF NEEDED)
#pretitle = "0 - "
pretitle = ""

fillvalue = "-999.0"

# Grab all the directory names from the CTD data directory.
ctd_dat = '/raid/baum/CTD_ERDDAP/CTD/'

errors = open(ctd_dat + "/000-PROBLEMS", "w")

#dirs = glob.glob(ctd_dat + '/*')


dirs = ['/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0001',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0002',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0003',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0046',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0064',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0065',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0066',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0084',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0090',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0093',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0096',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0099',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0108',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x134.073.0005',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x134.073.0014',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x134.073.0021',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.076.0003',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.076.0004',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.076.0009',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.076.0013',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.077.0002',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.077.0003',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.077.0005',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.078.0005',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.078.0025',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.078.0026',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.078.0027',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0008',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0009',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0010',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0011',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0012',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0013',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0014',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0015',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0017',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0018',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0048',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0049',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x138.079.0059',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0001',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0006',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0008',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0009',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0012',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x139.145.0013',
 '/raid/baum/CTD_ERDDAP/CTD/R2.x222.000.0011',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x260.000.0023',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x261.000.0001',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x261.000.0002',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x261.000.0004',
 '/raid/baum/CTD_ERDDAP/CTD/Y1.x030.000.0007',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0067',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x134.073.0024',
 '/raid/baum/CTD_ERDDAP/CTD/R1.x134.073.0043',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x260.000.0018',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x260.000.0047',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x260.204.0021',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x267.179.0019',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x267.179.0026',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x268.000.0017',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x268.000.0002',
 '/raid/baum/CTD_ERDDAP/CTD/R4.x268.000.0006']

# '/raid/baum/CTD_ERDDAP/CTD/R4.x268.000.0048']

# R1.x137.130.0012

#dirs = ['/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0001']
#	'/raid/baum/CTD_ERDDAP/CTD/R1.x132.134.0002']

# Iterate over all directory names.
for dirname in dirs:

	pathdirname = dirname
	dirname = os.path.basename(dirname)

	ncout = pathdirname + '.nc'

# ctd_prototype.xml
	erddap_template_xml = ctd_dat + "ERDDAP_PROFILE.XML"
	erddap_output_xml = ctd_dat + "erddap-" + dirname + ".xml"
#	print (" template = ",erddap_template_xml)
#	print (" output   = ",erddap_output_xml)

#  Scrape the title information from the GOMRI site.

	av,bv,cv,dv = dirname.split('.')
	udi = str(av + "." + bv + "." + cv + ":" + dv)

#	print (" Scraping: ",gomri_url+udi)
	res = requests.get(gomri_url+udi)
	res.raise_for_status()
	soup = bs4.BeautifulSoup(res.text,"lxml")
	title = "CTD - " + udi + " - " + soup.select('h2')[0].getText()
	title = pretitle + title.encode('utf-8')
#	print (" title = ",title)

#  Create ERDDAP XML section for this dataset from a template.

	replacements = {'DATASET_ID':dirname.replace('.','_'),
		        'NETCDF_FILENAME':os.path.basename(ncout),
			'SOURCE_URL':gomri_url+udi,
			'TITLE_STRING':title,
			'FILE_DIR':erddap_meta_filedir}

	with open(erddap_template_xml) as infile, open(erddap_output_xml, 'w') as outfile:
		for line in infile:
			for src, target in replacements.iteritems():
				line = line.replace(src, target)
			outfile.write(line)

#  Create THREDDS XML section for this dataset from a template.

	thredds_template_xml = ctd_dat + "THREDDS_PROFILE.XML"
	thredds_output_xml = ctd_dat + "thredds-" + dirname + ".xml"

	replacements = {'THREDDS_NAME':dirname,
			'THREDDS_ID':dirname,
			'THREDDS_PATH':"profile/"+dirname,
			'THREDDS_LOCATION':"/data/thredds/profile/"+dirname}

	with open(thredds_template_xml) as infile, open(thredds_output_xml, 'w') as outfile:
		for line in infile:
			for src, target in replacements.iteritems():
				line = line.replace(src, target)
			outfile.write(line)

	print (" Processing directory: ",dirname)
	errors.write ("Processing directory: " + dirname + "\n")

#  Create filename for output netCDF file and open it.

	filepath = ctd_dat + dirname + '/*.nc'
	files = glob.glob(filepath)

	ncout = pathdirname + '.nc'

	nc = netCDF4.Dataset(ncout,'w')

#  Create dimension(s) and variables for output netCDF file.

	nc.createDimension('station',len(files))
#        nc.createDimension('station',None)
        stat = nc.createVariable('station',str,'station')
        stat.long_name = "station identification number"

	nc.dataset_id = udi
	nc.griidc_udi = udi

        dap_url = nc.createVariable('dap_url',str,'station')
        dap_url.long_name = "THREDDS OpenDAP URL"
        http_url = nc.createVariable('http_url',str,'station')
        http_url.long_name = "THREDDS HTTP URL"

        otime = nc.createVariable('time','f8',('station',))
        otime.standard_name = "time"
        otime.long_name = "time"
        otime.calendar = "Julian"

        olon = nc.createVariable('lon','f8',('station',))
        olon.standard_name = "longitude"
        olon.long_name = "longitude"
        olon.units = "degrees_east"
        olon.axis = "X"

        olat = nc.createVariable('lat','f8',('station',))
        olat.standard_name = "latitude"
        olat.long_name = "latitude"
        olat.units = "degrees_north"
        olat.axis = "Y"

        odep = nc.createVariable('depth','f8',('station',),fill_value=fillvalue)
        odep.long_name = "depth"
        odep.standard_name = "depth"
        odep.units = "m"
        odep.positive = "down"
        odep.axis = "Z" 

	onum = nc.createVariable('vertical_points','i4',('station',))
	onum.long_name = "number of vertical points"
	onum.standard_name = "vertical_points"
	onum.units = "1"

#  Extract dataset-wide info from typical file, i.e. the first file.

	nc_in = netCDF4.Dataset(files[0], 'r')
	nc_dict = nc_in.__dict__

	consortia = nc_dict['project']
	nc.consortia = consortia
	nc.organization = "GRIIDC"
	nc.title = nc_dict['title']
	nc.summary = nc_dict['summary']
	try:
		nc.platform = nc_dict['platform']
	except:
		print (" No platform in: ", files[0])
	nc.instrument = nc_dict['instrument']
	nc.sea_name = nc_dict['sea_name']
	try:
		nc.originators_cruise_code = nc_dict['Originators_Cruise_Code']
	except:
		errstr = " No Originators_Cruise_Code in: " + files[0]
		print (errstr)
		errors.write (str(errstr) + "\n")
	try:
		nc.submitting_institution = nc_dict['submitting_institution']
	except:
		errstr = " No submitting_institution in: " + files[0]
		print (errstr)
		errors.write (str(errstr) + "\n")
	nc.publisher_name = nc_dict['publisher_name']
	nc.publisher_url = nc_dict['publisher_url']

#  Add list of variable names to meta-netCDF file.

	vert_variables = ",".join(nc_in.variables.keys())
	nc.dataset_variables = vert_variables

#  Add list of units to meta-netCDF file.

	try:
		string = ""
		for var in nc_in.variables.values():
			try:
				string = string + getattr(var,"units") + ","
			except:
				string = string + "NA" + ","
		nc.variable_units = string[:-1]
	except:
		errstr = " No unit attributes available in prototype file: " + os.path.basename(files[0])
		print (errstr)
		errors.write (errstr + "\n")

#  Add list of standard_name values ot meta-netCDF file.

	try:
		string = ""
		for var in nc_in.variables.values():
			try:
				string = string + getattr(var,"standard_name") + ","
			except:
				string = string + "NA" + ","
		nc.variable_standard_names = string[:-1]
	except:
		errstr = " No standard names available in prototype file: " + os.path.basename(files[0])
		print (errstr)
		errors.write (errstr + "\n")


#  Add list of long_name values to meta-netCDF file.

	try:
		string = ""
		for var in nc_in.variables.values():
			string = string + getattr(var,"long_name") + ","
		nc.variable_long_names = string[:-1]
	except:
		errstr = " No long names available in prototype file: " + os.path.basename(files[0])
		print (errstr)
		errors.write (errstr + "\n")

	nc_in.close()
	
#  Loop over all files in directory to snag per-file information.

	nf = 0
	for f in files:

#		errstr = " No long names available in prototype file: " + os.path.basename(files[0])print (" Opening file: ",f)

		try:
			nc_in = netCDF4.Dataset(f, 'r')
			nc_dict = nc_in.__dict__

		except:
			print (' Cannot open file: ',f)

		try:
			time = nc_in.variables["time"]
		except:
			print (" No time variable in: " + dirname + "/" + os.path.basename(f) + " or the time variable is not called 'time'.")
			print (" This dataset cannot be processed without a proper time variable.")
			del_string = "rm " + pathdirname + ".nc"
			print " del_string = ",del_string
			os.system(del_string)
			nc_in.close()
			errors.write(" ***** CRITICAL *****\n")
			errors.write(" No time variable in prototype file: " + os.path.basename(f) + ", or the time variable is not called 'time'.\n")
			errors.write(" This dataset cannot be processed without a proper time variable.\n")
			break	

		time_dict = time.__dict__
		time_units = time_dict["units"]
		otime.units = time_units

		time = nc_in.variables["time"][:]

		try:
			lat = nc_in.variables["lat"][:]
			olat[nf] = lat[0]
		except:
			print (" No latitude variable in: " + dirname + "/" + os.path.basename(f) + " or the latitude variable is not called 'lat'.")
			olat[nf] = fillvalue

		try:
			lon = nc_in.variables["lon"][:]
			olon[nf] = lon[0]
		except:
			print (" No longitude variable in: " + dirname + "/" + os.path.basename(f) + " or the longitude variable is not called 'lon'.")
			olat[nf] = fillvalue

		try:
			dep = nc_in.variables["dep"][:]
			depmax = max(dep)
			odep[nf] = depmax	
#			odep[nf] = dep[0]
		except:
			errstr = " No depth variable in prototype file: " + os.path.basename(f) + ", or the depth variable is not called 'dep'."
			print (errstr)
			errors.write (errstr + "\n")
			odep[nf] = fillvalue

		onum[nf] = nc_in.dimensions.values()[1].size
		otime[nf] = time[0]
		stat[nf] = os.path.basename(f)[:-3]

		dap = thredds_url + "/thredds/dodsC/profile/" + dirname + "/"
		http = thredds_url + "/thredds/fileServer/profile/" + dirname + "/"
		dap = dap + os.path.basename(f) + ".html"
		http = http + os.path.basename(f)

		http_url[nf] = http
		dap_url[nf] = dap

		nc_in.close()

		nf = nf + 1

	nc.close()
	print ("   ")
	errors.write("\n")

errors.close()

