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
            indices[i].checked = indices[i].checked != true;
    }
}


function checkNone(){
    var indices = document.getElementsByTagName("input");
    for (var i = 0; i < indices.length; i++){
        if(indices[i].type == "checkbox")
            indices[i].checked = false;
    }
}


function checkFileExtend(id){
    // Check the upload file type, it must be txt or fasta type.
    var filePath = document.getElementById(id).value;

    if (filePath != ""){
        var re = /(\\+)/g;
        filePath = filePath.replace(re,"#");
        var path_split = filePath.split("#");
        var filename = path_split[path_split.length - 1];
        var name_split = filename.split(".");
        var extend = name_split[name_split.length - 1];
        var extendAllowed = "txt, fasta";
        var resIndex = extendAllowed.lastIndexOf(extend);
        if (resIndex >= 0)
            return true;
        else{
            alert("The upload file must be txt or fasta type.");
            if (id == "upload_data")
                document.myForm.upload_data.focus();
            else if (id == "upload_ind")
                document.myForm.upload_ind.focus();

            return false;
        }
    }

    return true;
}


function formCheckArgs(mode)
{
    // Check parameter lag.
    if (mode == "Auto covariance" || mode == "Cross covariance" || mode == "Auto-cross covariance")
        if(document.myForm.lag.value == ""){
            alert("Please input the parameter lag!");
            document.myForm.lag.focus();
            return false;
        }else if(document.myForm.lag.value <= 0){
            alert("Parameter lag must be larger than 0!");
            document.myForm.lag.focus();
            return false;
        }

    // Check parameter λ, w.
    if (mode == "PseSSC" || mode == "PseDPC" || mode == "pPseDNC" || mode == "sPseDNC"){
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

    // Check upload file format.
    if (checkFileExtend("upload_data") == false)
        return false;
    try{
        if (checkFileExtend("upload_ind") == false)
            return false;
    } catch (err){
    }

    // Physicochemical properties and index file cannot be both null.
    if (mode == "Auto covariance" || mode == "Cross covariance" || mode == "Auto-cross covariance" || mode == "pPseDNC" || mode == "sPseDNC"){
        if (document.getElementById("upload_ind").value == ""){
            var checked = false;
            var indices = document.getElementsByTagName("input");
            for (var i = 0; i < indices.length; i++){
                if(indices[i].type == "checkbox" && indices[i].checked == true){
                    checked = true;
                    break;
                }
            }
            if (checked == false){
                alert("The physicochemical properties checkbox and user-defined physicochemical properties index file cannot be both null!");
                document.myForm.upload_ind.focus();
                return false;
            }
        }
    }

    // Sequence input cannot be both null.
    if (document.getElementById("rec_data").value == "" && document.getElementById("upload_data").value == ""){
        alert("You must input the sequences in textarea or upload sequence file.");
        document.myForm.rec_data.focus();
        return false;
    }

    // Sequence input cannot be both true.
    if (document.getElementById("rec_data").value != "" && document.getElementById("upload_data").value != ""){
        alert("You cannot both input the sequences in textarea and upload sequence file.");
        document.myForm.rec_data.focus();
        return false;
    }

    return true;
}


function setAutoForm(mode){
    if(mode == "Auto covariance" || mode == "Cross covariance" || mode == "Auto-cross covariance" || mode == "pPseDNC" || mode == "sPseDNC"){
        var indices = document.getElementsByTagName("input");
        for (var i = 0; i < indices.length; i++){
            if (indices[i].type == "checkbox"){
                indices[i].checked = !!(indices[i].name == "Roll (RNA)" || indices[i].name == "Rise (RNA)"
                || indices[i].name == "Shift (RNA)" || indices[i].name == "Slide (RNA)"
                || indices[i].name == "Tilt (RNA)" || indices[i].name == "Twist (RNA)");
            }
        }

        if (mode == "Auto covariance" || mode == "Cross covariance" || mode == "Auto-cross covariance"){
            document.forms[0].lag.value = 10;
        }else{
            document.forms[0].lamada.value = 10;
            document.forms[0].w.value = 0.05;
        }
    }else if (mode == "PseSSC"){
        document.forms[0].n.value = 2;
        document.forms[0].lamada.value = 13;
        document.forms[0].w.value = 0.5;
    }else if (mode == "PseDPC"){
        document.forms[0].n.value = 7;
        document.forms[0].lamada.value = 15;
        document.forms[0].w.value = 1;
    }

    document.forms[0].rec_data.value = ">EXAMPLE1\nGCAUCCGGGUUGAGGUAGUAGGUUGUAUGGUUUAGAGUUACACCCUGGGAGUUAACUGUACAACCUUCUAGCUUUCCUUGGAGC\n>EXAMPLE2\nCCUAGGAAGAGGUAGUAGGUUGCAUAGUUUUAGGGCAGGGAUUUUGCCCACAAGGAGGUAACUAUACGACCUGCUGCCUUUCUUAGG\n>EXAMPLE3\nGAGGGCAGGGGGCACAGUCCAACUCCAGGCUUGUAGCUGUCCAGGGGCUGGGUGCCCGCCCGGCAGCGGCAGACUGUGUCCUGUGUGGCCGUGCACA\n>EXAMPLE4\nAAGCACCAGAGACGAACAGUCUGGUGUCAGGAGCAAGAGAAAGGCUCAUGACAUCUCCAGUGUGUCCGGUAAACGUGGUCG";
}