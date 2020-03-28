
# create a configuration object that changes the preprocessors
from nbconvert import HTMLExporter
from traitlets.config import Config

c = get_config()

c.NbConvertApp.export_format = 'html'

c.NbConvertApp.output_files_dir = './html'
#c.NbConvertApp.output_files_dir = './html/FCI_Demo_KnockdownAndGrind2'

c.HTMLExporter.preprocessors = [
	'nbconvert.preprocessors.ExecutePreprocessor',
	'nbconvert.preprocessors.coalesce_streams',
	'nbconvert.preprocessors.ExtractOutputPreprocessor']

c.ExecutePreprocessor.allow_errors = True

#c.FilesWriter.build_directory = './html/FCI_Demo_StandConversion_AlderToConifer'
c.FilesWriter.build_directory = './html'





