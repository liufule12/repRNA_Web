/**
 * Created by aleeee on 15-1-19.
 */


function checkAll(){
    var indices = document.getElementsByTagName("input");
    for (var i = 0; i < indices.length; i++) {
        if(indices[i].type == "checkbox")
            indices[i].checked = true;
    }
}


function checkRev(){
    var indices = document.getElementsByTagName("input");
    for (var i = 0; i < indices.length; i++){
        if(indices[i].type == "checkbox")
            if(indices[i].checked == true)
                indices[i].checked = false;
            else
                indices[i].checked = true;
    }
}


function checkNone(){
    var indices = document.getElementsByTagName("input");
    for (var i = 0; i < indices.length; i++){
        if(indices[i].type == "checkbox")
            indices[i].checked = false;
    }
}


//function checkFileExtend(filename){
//    var extend = filename.split(".")[-1];
//    alert(extend);
//}


function formCheckArgs(mode)
{
    if (mode == "DAC" || mode == "DCC" || mode == "DACC")
        if(document.myForm.lag.value == ""){
            alert("Please input the parameter lag!");
            document.myForm.lag.focus();
            return false;
        }else if(document.myForm.lag.value < 0){
            alert("Parameter lag must be larger than 0!");
            document.myForm.lag.focus();
            return false;
        }

    if (mode == "PseSSC" || mode == "PseDPC" || mode == "PC-PseDNC-General" || mode == "SC-PseDNC-General"){
        if (document.myForm.lamada.value == ""){
            alert("Please input the parameter λ!");
            document.myForm.lamada.focus();
            return false;
        }else if(document.myForm.lamada.value < 0){
            alert("Parameter λ must be larger than 0!");
            document.myForm.lamada.focus();
            return false;
        }

        if (document.myForm.w.value == ""){
            alert("Please input the parameter w!");
            document.myForm.w.focus();
            return false;
        }else{
            try{
                var w = parseFloat(document.myForm.w.value);
                if (w < 0 || w > 1){
                    alert("Parameter w must be ranging from 0 to 1!");
                    document.myForm.w.focus();
                    return false;
                }
            }
            catch (err){
                alert("Parameter w must be range from 0 to 1!");
                document.myForm.w.focus();
                return false;
            }
        }
    }

    //var filename = fileExtend(document.myForm.upload_data.value);
    //alert(filename);
    //var extend = filename.split(".")[-1];
    //alert(extend);

    return(true);
}