
def getFieldValues(shapefile,target_field):
    ##returns the number of values in an attribute table
    import arcpy, sys, os

    rows=arcpy.SearchCursor(shapefile,"","",target_field)
    value_list =[]
    for row in rows:
        this_value = row.getValue(target_field)
        if this_value not in value_list:
            value_list.append(this_value)        
    return value_list
        
