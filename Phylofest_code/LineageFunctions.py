
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
        
def getFieldMinMax(layer,target_field):
    ##returns the minimum and maximum values in a numeric field of an attribute table
    import arcpy.analysis

    #fields = [target_field]
    rows=arcpy.da.SearchCursor(layer,(target_field))
    #value_list =[row[0] for row in rows]
    
    value_list=[]
    for row in rows:
        #this_value = row.getValue(target_field)
        this_value = row[0]
        if this_value not in value_list:
            if this_value != None:
                value_list.append(this_value)
    minval=min(value_list)
    maxval=max(value_list)
    return [minval,maxval]
