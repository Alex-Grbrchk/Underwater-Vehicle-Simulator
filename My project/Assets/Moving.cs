using System.Collections;
using System.Collections.Generic;
using UnityEngine;
using UnityEngine.UI;

public class Moving : MonoBehaviour
{
    public float height, speed, HForce = 2f;
    public Rigidbody RB;

    // Start is called before the first frame update
    void Start()
    {
        speed = 0f;
        height = 0f;
        RB = gameObject.GetComponent<Rigidbody>();
    }

    // Update is called once per frame
    private void FixedUpdate()
    {
        if (gameObject.transform.position.y != height)
        {
            RB.AddForce(new Vector3(0, (height-gameObject.transform.position.y)*HForce * Time.fixedDeltaTime, 0));
            if(RB.velocity.y > Mathf.Abs(height - gameObject.transform.position.y))
            {
                RB.velocity = new Vector3(RB.velocity.x, Mathf.Abs(height - gameObject.transform.position.y), RB.velocity.z);
            }
        }


    }
    void Update()
    {

        if (Input.GetKey(KeyCode.UpArrow))
        {
            this.transform.Translate(new Vector3(0, 0, speed) * Time.deltaTime);
        }
        if (Input.GetKey(KeyCode.DownArrow)) 
        {
            this.transform.Translate(new Vector3(0, 0, -speed) * Time.deltaTime);
        }
        if (Input.GetKey(KeyCode.LeftArrow))
        {
            this.transform.Rotate(0, -30 * Time.deltaTime, 0);

        }
        if (Input.GetKey(KeyCode.RightArrow))
        {
            this.transform.Rotate(0, 30 * Time.deltaTime, 0);
        }
    }
}
