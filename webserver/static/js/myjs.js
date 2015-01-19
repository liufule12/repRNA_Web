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

