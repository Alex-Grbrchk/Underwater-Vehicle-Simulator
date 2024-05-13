using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class ViewChage : MonoBehaviour
{
    // Start is called before the first frame update\
    public GameObject Cam1, Cam2, Cam1s, Cam2s;
    public Button View;
    public bool f = true;
    
    void Awake()
    {
        Cam1s.SetActive (true);
        Cam2.SetActive (true);
        Cam1.SetActive(false);
        Cam2s.SetActive(false);
        View.onClick.AddListener(ViewCh);
 
    }

    // Update is called once per frame
    void Update()
    {
        if (Input.GetKeyDown(KeyCode.Q))
        {
            ViewCh();
          
        }
    }

  

    void ViewCh()
    {
        if (Cam1s.activeSelf == true)
        {
            Cam1s.SetActive(false);
            Cam2.SetActive(false);
            Cam1.SetActive(true);
            Cam2s.SetActive(true);
        }
        else
        {
            Cam1s.SetActive(true);
            Cam2.SetActive(true);
            Cam1.SetActive(false);
            Cam2s.SetActive(false);
        }
    }
   
}
