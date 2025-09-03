'''
Function that wraps the necessary methods from morgana to analyse images with
one or multiple gastruloids.
'''
import os
import sys

from morgana.ImageTools.objectsparsing import objectsparser
from morgana.DatasetTools.morphology import computemorphology
from morgana.DatasetTools.fluorescence import computefluorescence
from morgana.DatasetTools.fluorescence import io as ioFluo
from morgana.DatasetTools.morphology import io as ioMorph
from morgana.DatasetTools.straightmorphology import computestraightmorphology
from morgana.DatasetTools.straightmorphology import io as ioStraightMorph

# image_folder = r'Y:\Users\Apolline\time_lapse_experiments\gastruloids\expg83\d04\results\05_properties\images'
# mask_folder = r'Y:\Users\Apolline\time_lapse_experiments\gastruloids\expg83\d04\results\05_properties\masks'
# identifier_string = ''  
# objects_at_border = False
# objectsparser.parsing_images_touching_labels(image_folder, mask_folder, identifier_string, objects_at_border) # THIS IS WHAT GUILLERMO USES
# # extract data from all the folders
# input_folder = os.path.join(image_folder, 'splitObjects')
# save_folder = os.path.join(input_folder,'result_segmentation')
# data_fluo = computefluorescence.compute_fluorescence_info(input_folder)

def run_morgana_pipeline(image_folder, mask_folder, cond):
    identifier_string = ''
    objects_at_border = True

    print(image_folder)
    print(mask_folder)

    objectsparser.parsing_images_touching_labels(image_folder, mask_folder, identifier_string, objects_at_border) # THIS IS WHAT GUILLERMO USES
    #objectsparser.parsing_images(image_folder, mask_folder, identifier_string, objects_at_border) # THIS SHOULD WORK FOR MARTA


    # extract data from all the folders
    input_folder = os.path.join(image_folder, 'splitObjects')
    save_folder = os.path.join(input_folder,'result_segmentation')

    #maskType = "Unprocessed" #Straightened
    isTimelapse=False

    #if maskType == "Unprocessed":
    #data = computemorphology.compute_morphological_info(input_folder)
    #save_morphological_info = ioMorph.save_morpho_params( save_folder, cond, data )

    data_fluo = computefluorescence.compute_fluorescence_info(input_folder)
    save_fluo_info = ioFluo.save_fluo_info(save_folder, cond, data_fluo)
    
    #else:
    #data = computestraightmorphology.compute_straightwor_morphological_info(input_folder)
    #save_morphological_info = ioStraightMorph.save_straight_morpho_params( save_folder, cond, data )

if __name__ == '__main__':
    run_morgana_pipeline(sys.argv[1], sys.argv[2], sys.argv[3])

    # example
    # python \\larch\gs714\Desktop\run_morgana_pipeline.py Y:\Users\Guillermo\experiments_imaging\gastruloids\exp0013\results\test_morgana\images Y:\Users\Guillermo\experiments_imaging\gastruloids\exp0013\results\test_morgana\masks image_name
